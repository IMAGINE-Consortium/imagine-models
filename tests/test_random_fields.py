import ImagineModels as img

import pytest
import numpy as np


pytestmark = pytest.mark.skipif(not img.__has_random_fields__, reason="ImagineModels was built without FFTW")

vector_models = ['JF12RandomField', 'ESRandomField']
scalar_models = ['GaussianScalarField', 'LogNormalScalarField']

grid = img.RegularGrid(shape=[16, 16, 16], reference_point=[-1., -1., -1.], increment=[.125, .125, .125])


def _sample(model_string, seed):
    return getattr(img, model_string)().sample(grid, seed)


@pytest.mark.parametrize('model_string', vector_models + scalar_models)
def test_sample_smoke(model_string):
    field = _sample(model_string, seed=3)
    expected_shape = (3, *grid.shape) if model_string in vector_models else tuple(grid.shape)
    assert field.shape == expected_shape
    assert np.isfinite(field).all()
    assert np.std(field) > 0.


@pytest.mark.parametrize('model_string', vector_models + scalar_models)
def test_seed_reproducibility(model_string):
    assert np.array_equal(_sample(model_string, seed=7), _sample(model_string, seed=7))
    assert not np.array_equal(_sample(model_string, seed=7), _sample(model_string, seed=8))


@pytest.mark.parametrize('model_string', vector_models + scalar_models)
def test_profile_and_random_numbers(model_string):
    model = getattr(img, model_string)()
    irregular = img.IrregularGrid(np.asarray([0., 1.]), np.asarray([2.]), np.asarray([-1., 0., 1.]))
    profile = model.profile(irregular)
    assert profile.shape == (2, 1, 3)
    numbers = model.random_numbers(grid, 3)
    assert numbers.shape == ((3, *grid.shape) if model_string in vector_models else tuple(grid.shape))


def test_sample_requires_regular_grid():
    irregular = img.IrregularGrid(np.asarray([0., 1.]), np.asarray([2.]), np.asarray([-1., 0., 1.]))
    with pytest.raises(TypeError):
        img.GaussianScalarField().sample(irregular, 3)


def test_profile_matches_spatial_profile():
    model = img.JF12RandomField()
    irregular = img.IrregularGrid(np.asarray([-8.5, 3.]), np.asarray([0., 4.]), np.asarray([0.1]))
    profile = model.profile(irregular)
    assert profile[1, 1, 0] == model.spatial_profile(3., 4., .1)
