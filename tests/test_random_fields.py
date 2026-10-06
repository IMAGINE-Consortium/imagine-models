import ImagineModels as img

import pytest
import numpy as np


pytestmark = pytest.mark.skipif(not img.__has_random_fields__, reason="ImagineModels was built without FFTW")

vector_models = ['JF12RandomField', 'ESRandomField']
scalar_models = ['GaussianScalarField', 'LogNormalScalarField']

shape = [16, 16, 16]
increment = [.125, .125, .125]
zeropoint = [-1., -1., -1.]


def _evaluate(model_string, seed):
    mo = getattr(img, model_string)()
    return np.asarray(mo.on_grid(shape=shape, reference_point=zeropoint, increment=increment, seed=seed))


@pytest.mark.parametrize('model_string', vector_models + scalar_models)
def test_on_grid_smoke(model_string):
    field = _evaluate(model_string, seed=3)
    expected_shape = (3, *shape) if model_string in vector_models else tuple(shape)
    assert field.shape == expected_shape
    assert np.isfinite(field).all()
    assert np.std(field) > 0.


@pytest.mark.parametrize('model_string', vector_models + scalar_models)
def test_seed_reproducibility(model_string):
    assert np.array_equal(_evaluate(model_string, seed=7), _evaluate(model_string, seed=7))
    assert not np.array_equal(_evaluate(model_string, seed=7), _evaluate(model_string, seed=8))
