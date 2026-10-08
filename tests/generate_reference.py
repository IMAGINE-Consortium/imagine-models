import argparse
import json
import subprocess
import sys

import numpy as np

import reference_cases as rc


def _commit():
    return subprocess.run(["git", "rev-parse", "HEAD"], cwd=rc.REPO, capture_output=True, text=True).stdout.strip()


def _write(case, arrays, meta):
    rc.DATA_DIR.mkdir(exist_ok=True)
    np.savez_compressed(rc.DATA_DIR / f"{case}.npz", meta=json.dumps(meta), **arrays)


def regular(case, name, settings, commit):
    model = rc.make_model(name, settings)
    arrays = {"positions": rc.POSITIONS, "at_position": rc.at_positions(model, rc.POSITIONS)}
    if case not in rc.NO_GRID_CASES:
        arrays["on_regular_grid"] = rc.on_regular_grid(model)
        arrays["on_irregular_grid"] = rc.on_irregular_grid(model)
    meta = {
        "case": case,
        "model": name,
        "settings": settings,
        "commit": commit,
        "defaults": rc.parameter_defaults(model),
    }
    if rc.has_jacobian(model):
        columns = rc.jacobian_column_labels(name, model)
        arrays["jacobian_positions"] = rc.JACOBIAN_POSITIONS
        arrays["jacobian"] = np.array([rc.jacobian(model, p, columns) for p in rc.JACOBIAN_POSITIONS])
        meta["jacobian_columns"] = columns
        meta["jacobian_known_bad"] = rc.KNOWN_BAD_JACOBIAN.get(name)
    _write(case, arrays, meta)


def random(case, name, settings, commit):
    arrays = {}
    for grid_name, grid in rc.RANDOM_GRIDS.items():
        for seed in rc.RANDOM_SEEDS:
            arrays[f"sample_{grid_name}_{seed}"] = rc.sample(rc.make_model(name, settings), grid, seed)
    meta = {
        "case": case,
        "model": name,
        "settings": settings,
        "commit": commit,
        "defaults": rc.parameter_defaults(rc.make_model(name, settings)),
    }
    _write(case, arrays, meta)


def main():
    parser = argparse.ArgumentParser(description="Generate reference data for tests/test_reference.py")
    parser.add_argument("--force", action="store_true", help="overwrite existing reference data")
    parser.add_argument("cases", nargs="*", help="only (re)generate these cases")
    args = parser.parse_args()
    if not rc.RANDOM_CASES:
        sys.exit("ImagineModels was built without FFTW; cannot generate random-field reference data")
    selected = args.cases or list(rc.REGULAR_CASES) + list(rc.RANDOM_CASES)
    unknown = [c for c in selected if c not in rc.REGULAR_CASES and c not in rc.RANDOM_CASES]
    if unknown:
        sys.exit(f"unknown cases: {unknown}")
    existing = [c for c in selected if (rc.DATA_DIR / f"{c}.npz").exists()]
    if existing and not args.force:
        sys.exit(f"{len(existing)} cases already exist; reference data must not be regenerated casually (use --force)")
    commit = _commit()
    for case in selected:
        if case in rc.REGULAR_CASES:
            regular(case, *rc.REGULAR_CASES[case], commit)
        else:
            random(case, *rc.RANDOM_CASES[case], commit)
    print(f"wrote {len(selected)} cases to {rc.DATA_DIR}")


if __name__ == "__main__":
    main()
