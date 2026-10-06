import gc

import ImagineModels as img

import pytest
import numpy as np


regular_models = ['JaffeMagneticField', 'HelixMagneticField', 'JF12RegularField', 'SunMagneticField', 'UFMagneticField']
scalar_models = ['YMW16', 'UniformDensityField']

known_positions = {'JF12RegularField': [([0, 0, 0], [0, 0, 0]), # Galactic center
                                        ([.1, .3, .4], [0, 0, 0]), # Within inner boundary
                                        ([3, 20., .3], [0, 0, 0, ]),], # outside outer boundary
                   'JaffeMagneticField': [([0, 0, 0], [0., 0., 0.]),], # origin
                   'HelixMagneticField': [([0, 0, 0], [0, 0, 0]), # origin
                                         ],
                   'SunMagneticField': [([0, 0, 0], [0.41582338163551863, 1.9562952014676114, 0.0]), # origin
                                         ],
                   'UFMagneticField': [],
                   }

regular_grid = img.RegularGrid(shape=[2, 3, 4], reference_point=[-2., 3., .1], increment=[.1, 3., .1])
irregular_grid = img.IrregularGrid(np.asarray([-2.2, -1., 0., .1]), np.asarray([-1., 0., 3.]), np.asarray([-5., 0., .1]))


def _as_irregular(grid):
    axes = [grid.reference_point[d] + np.arange(grid.shape[d]) * grid.increment[d] for d in range(3)]
    return img.IrregularGrid(*axes)


def test_at_position():
    for model_string in regular_models:
        mo = getattr(img, model_string)()
        for position, value in known_positions[model_string]:
            pos = mo.at_position(*position)
            for j in range(3):
                assert pos[j] == value[j]


@pytest.mark.parametrize('model_string', regular_models + scalar_models)
def test_grid_consistency(model_string):
    mo = getattr(img, model_string)()
    on_regular = mo.evaluate(regular_grid)
    assert np.array_equal(on_regular, mo.evaluate(_as_irregular(regular_grid)), equal_nan=True)

    on_irregular = mo.evaluate(irregular_grid)
    for i, x in enumerate(irregular_grid.x):
        for j, y in enumerate(irregular_grid.y):
            for k, z in enumerate(irregular_grid.z):
                expected = np.asarray(mo.at_position(x, y, z), dtype=float)
                assert np.array_equal(on_irregular[..., i, j, k], expected, equal_nan=True)


def test_return_shapes():
    assert img.JF12RegularField().evaluate(regular_grid).shape == (3, 2, 3, 4)
    assert img.JF12RegularField().evaluate(irregular_grid).shape == (3, 4, 3, 3)
    assert img.YMW16().evaluate(regular_grid).shape == (2, 3, 4)
    assert img.YMW16().evaluate(irregular_grid).shape == (4, 3, 3)


def test_returned_array_owns_its_memory():
    mo = img.UniformMagneticField()
    mo.bx = 1.5
    field = mo.evaluate(regular_grid)
    del mo
    gc.collect()
    assert np.all(field[0] == 1.5)
    assert np.all(field[1:] == 0.)


def test_grid_interface():
    assert regular_grid.shape == [2, 3, 4]
    assert regular_grid.size == 24
    assert irregular_grid.shape == [4, 3, 3]
    assert np.array_equal(irregular_grid.y, [-1., 0., 3.])
    img.RegularGrid(np.asarray([2, 3, 4]), [-2., 3., .1], [.1, 3., .1])

    with pytest.raises(TypeError):
        img.RegularGrid([2., 3., 4.], [-2., 3., .1], [.1, 3., .1]) # shape must be int
    with pytest.raises(TypeError):
        img.RegularGrid([2, 3, 4.], [-2., 3., .1], [.1, 3., .1]) # all shape must be int
    with pytest.raises(ValueError):
        img.RegularGrid([2, 0, 4], [-2., 3., .1], [.1, 3., .1]) # shape must be positive
    with pytest.raises(ValueError):
        img.IrregularGrid(np.asarray([1., 2.]), np.asarray([]), np.asarray([0.])) # axes must not be empty
    with pytest.raises(TypeError):
        img.UniformMagneticField().evaluate([2, 3, 4]) # a grid object is required


def test_parameter_update():
    umf = img.UniformMagneticField()
    assert umf.bx == 0.
    assert umf.by == 0.
    assert umf.bz == 0.

    assert umf.at_position(2.4, 2.1, -.2) == (0., 0., 0.)

    umf.bx = -3.2

    assert umf.bx == -3.2

    assert umf.at_position(2.4, 2.1, -.2) == (-3.2, 0., 0.)
