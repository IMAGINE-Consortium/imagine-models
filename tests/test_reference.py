import json

import numpy as np
import pytest

import reference_cases as rc


def _load(case):
    path = rc.DATA_DIR / f"{case}.npz"
    if not path.exists():
        pytest.skip(f"no reference data for {case}")
    with np.load(path) as f:
        data = {key: f[key] for key in f.files}
    return data, json.loads(str(data["meta"]))


def _assert_close(actual, expected, rtol, atol):
    assert actual.shape == expected.shape
    np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol, equal_nan=True)


def _has_key(case, key):
    path = rc.DATA_DIR / f"{case}.npz"
    if not path.exists():
        return False
    with np.load(path) as f:
        return key in f.files


def _cases_with(key):
    return [c for c in rc.REGULAR_CASES if _has_key(c, key)]


@pytest.mark.parametrize("case", rc.REGULAR_CASES)
def test_at_position(case):
    data, meta = _load(case)
    model = rc.make_model(meta["model"], meta["settings"])
    _assert_close(rc.at_positions(model, data["positions"]), data["at_position"], rc.RTOL_REGULAR, rc.ATOL_REGULAR)


@pytest.mark.parametrize("case", _cases_with("on_regular_grid"))
def test_on_regular_grid(case):
    data, meta = _load(case)
    model = rc.make_model(meta["model"], meta["settings"])
    _assert_close(rc.on_regular_grid(model), data["on_regular_grid"], rc.RTOL_REGULAR, rc.ATOL_REGULAR)


@pytest.mark.parametrize("case", _cases_with("on_irregular_grid"))
def test_on_irregular_grid(case):
    data, meta = _load(case)
    model = rc.make_model(meta["model"], meta["settings"])
    _assert_close(rc.on_irregular_grid(model), data["on_irregular_grid"], rc.RTOL_REGULAR, rc.ATOL_REGULAR)


@pytest.mark.parametrize("case", rc.REGULAR_CASES)
def test_parameter_defaults(case):
    _, meta = _load(case)
    model = rc.make_model(meta["model"], meta["settings"])
    assert rc.parameter_defaults(model) == meta["defaults"]


@pytest.mark.parametrize("case", _cases_with("jacobian"))
def test_jacobian(case):
    data, meta = _load(case)
    model = rc.make_model(meta["model"], meta["settings"])
    if not rc.has_jacobian(model):
        pytest.skip("ImagineModels was built without autodiff")
    actual = np.array([rc.jacobian(model, p, meta["jacobian_columns"]) for p in data["jacobian_positions"]])
    _assert_close(actual, data["jacobian"], rc.RTOL_REGULAR, rc.ATOL_REGULAR)


def _finite_difference(model, label, position, h):
    p0 = getattr(model, label)
    setattr(model, label, p0 + h)
    up = np.atleast_1d(np.asarray(model.at_position(*position), dtype=float))
    setattr(model, label, p0 - h)
    down = np.atleast_1d(np.asarray(model.at_position(*position), dtype=float))
    setattr(model, label, p0)
    return (up - down) / (2 * h)


def _one_sided_difference(model, label, position, h):
    p0 = getattr(model, label)
    center = np.atleast_1d(np.asarray(model.at_position(*position), dtype=float))
    setattr(model, label, p0 + h)
    shifted = np.atleast_1d(np.asarray(model.at_position(*position), dtype=float))
    setattr(model, label, p0)
    return (shifted - center) / h


def _column(jac, index, n_columns):
    jac = jac.reshape(-1, n_columns) if jac.size % n_columns == 0 else jac
    return jac[:, index]


def _finite_difference_cases():
    return [
        pytest.param(c, marks=pytest.mark.xfail(strict=True, reason=f"issue {rc.BROKEN_DERIVATIVE_CASES[c]}"))
        if c in rc.BROKEN_DERIVATIVE_CASES
        else c
        for c in _cases_with("jacobian")
    ]


@pytest.mark.parametrize("case", _finite_difference_cases())
def test_jacobian_matches_finite_differences(case):
    data, meta = _load(case)
    model = rc.make_model(meta["model"], meta["settings"])
    labels = meta["jacobian_columns"]
    failures = []
    for position, jac in zip(data["jacobian_positions"], data["jacobian"]):
        for index, label in enumerate(labels):
            scale = max(1.0, abs(getattr(model, label)))
            fine = _finite_difference(model, label, position, 1e-6 * scale)
            coarse = _finite_difference(model, label, position, 1e-4 * scale)
            if not np.allclose(fine, coarse, rtol=1e-2, atol=1e-6, equal_nan=True):
                continue
            expected = _column(jac, index, len(labels))
            atol = 1e-6 * max(1.0, np.nanmax(np.abs(fine), initial=0.0))
            candidates = [fine] + [_one_sided_difference(model, label, position, s * 1e-6 * scale) for s in (1, -1)]
            if not any(np.allclose(expected, c, rtol=1e-4, atol=atol, equal_nan=True) for c in candidates):
                failures.append((label, position.round(2).tolist(), expected.tolist(), fine.tolist()))
    assert not failures, failures


@pytest.mark.parametrize("case", rc.RANDOM_CASES)
@pytest.mark.parametrize("grid_name", rc.RANDOM_GRIDS)
@pytest.mark.parametrize("seed", rc.RANDOM_SEEDS)
def test_random_sample(case, grid_name, seed):
    data, meta = _load(case)
    model = rc.make_model(meta["model"], meta["settings"])
    actual = rc.sample(model, rc.RANDOM_GRIDS[grid_name], seed)
    _assert_close(actual, data[f"sample_{grid_name}_{seed}"], rc.RTOL_RANDOM, rc.ATOL_RANDOM)
