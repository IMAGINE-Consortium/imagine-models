import json
from pathlib import Path

import numpy as np
import pytest

import ImagineModels as img
import reference_cases as rc

DATA_DIR = Path(__file__).resolve().parent / "external_data"
CASES = sorted(p.stem for p in DATA_DIR.glob("*.npz"))


def _load(case):
    with np.load(DATA_DIR / f"{case}.npz") as f:
        data = {key: f[key] for key in f.files}
    return data, json.loads(str(data["meta"]))


def _evaluate(model, quantity, positions):
    if quantity == "rms":
        return np.asarray(model.rms(positions[:, 0], positions[:, 1], positions[:, 2]), dtype=float)[:, None]
    return rc.at_positions(model, positions)


def test_cases_present():
    assert len(CASES) > 0


@pytest.mark.parametrize("case", CASES)
def test_external_reference(case):
    data, meta = _load(case)
    if not hasattr(img, meta["model"]):
        pytest.skip(f"{meta['model']} not available in this build")
    assert np.count_nonzero(data["values"]) > data["values"].size // 4
    model = rc.make_model(meta["model"], meta["settings"])
    actual = _evaluate(model, meta["quantity"], data["positions"])
    np.testing.assert_allclose(actual, data["values"], rtol=meta["rtol"], atol=meta["atol"],
                               err_msg=f"{case}: {meta['source']}")
