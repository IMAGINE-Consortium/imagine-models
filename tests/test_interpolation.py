import numpy as np
import pytest

import ImagineModels as img

GRID = img.RegularGrid(shape=[5, 4, 3], reference_point=[-2.0, 1.0, -0.5], increment=[0.5, 0.25, 0.4])


def _axes(grid):
    return [grid.reference_point[a] + grid.increment[a] * np.arange(grid.shape[a]) for a in range(3)]


def _affine(x, y, z):
    return np.array([1.0 + 2.0 * x - 3.0 * y + 0.5 * z, x + y, 0.3 * z - x])


def _nodes(grid):
    return [c.ravel() for c in np.meshgrid(*_axes(grid), indexing="ij")]


def test_linear_reproduces_affine_fields():
    data = _affine(*np.meshgrid(*_axes(GRID), indexing="ij"))
    rng = np.random.default_rng(1)
    x, y, z = (
        rng.uniform(lo, lo + inc * (n - 1), 50) for lo, inc, n in zip(GRID.reference_point, GRID.increment, GRID.shape)
    )
    np.testing.assert_allclose(img.interpolate(data, GRID, x, y, z), _affine(x, y, z), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(img.interpolate(data[0], GRID, x, y, z), _affine(x, y, z)[0], rtol=1e-12, atol=1e-12)


def test_point_cloud_and_broadcasting_agree():
    data = img.JF12MagneticField().evaluate(GRID)
    x = np.linspace(-2.0, 0.0, 7)
    cloud = img.PointCloud(x, np.full_like(x, 1.3), np.zeros_like(x))
    np.testing.assert_array_equal(img.interpolate(data, GRID, cloud), img.interpolate(data, GRID, x, 1.3, 0.0))
    assert img.interpolate(data, GRID, np.zeros((2, 3)) - 1.0, 1.2, 0.0).shape == (3, 2, 3)
    assert img.interpolate(data[1], GRID, np.zeros((2, 3)) - 1.0, 1.2, 0.0).shape == (2, 3)


@pytest.mark.parametrize("method", ["linear", "nearest"])
def test_grid_nodes_return_grid_values(method):
    data = img.JF12MagneticField().evaluate(GRID)
    values = img.interpolate(data, GRID, *_nodes(GRID), method=method)
    np.testing.assert_allclose(values, data.reshape(3, -1), rtol=1e-12, atol=1e-14)


def test_random_sample_at_grid_nodes():
    if not img.__has_random_fields__:
        pytest.skip("ImagineModels was built without FFTW")
    grid = img.RegularGrid(shape=[8, 8, 4], reference_point=[-1.0, -1.0, -0.5], increment=[0.25, 0.25, 0.25])
    sample = img.UF26RandomField().sample(grid, 3)
    np.testing.assert_allclose(img.interpolate(sample, grid, *_nodes(grid)), sample.reshape(3, -1), rtol=1e-12)


def test_outside_and_invalid_input():
    data = img.JF12MagneticField().evaluate(GRID)
    with pytest.raises(img.GridError):
        img.interpolate(data, GRID, 0.5, 1.2, 0.0)
    assert np.all(np.isnan(img.interpolate(data, GRID, [0.5, -1.0], 1.2, [0.0, np.nan], nan_outside=True)))
    with pytest.raises(img.GridError):
        img.interpolate(data, GRID, -1.0, 1.2, 0.0, method="cubic")
    with pytest.raises(img.GridError):
        img.interpolate(data[:, :-1], GRID, -1.0, 1.2, 0.0)
    with pytest.raises(img.GridError):
        img.interpolate(data[0, 0], GRID, -1.0, 1.2, 0.0)
