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


def test_python_subclasses():
    class Linear(img.RegularScalarField):
        def at_position(self, x, y, z):
            return x + 2 * y + 3 * z

    class Constant(img.RegularVectorField):
        def at_position(self, x, y, z):
            return [1., 2., 3.]

    class Missing(img.RegularVectorField):
        pass

    assert np.array_equal(Linear().evaluate(img.IrregularGrid([1., 2.], [0.], [1.])).ravel(), [4., 5.])
    field = Constant().evaluate(regular_grid)
    assert field.shape == (3, 2, 3, 4)
    assert np.all(field[2] == 3.)
    with pytest.raises(RuntimeError):
        Missing().evaluate(regular_grid)


@pytest.mark.parametrize('model_string', regular_models + scalar_models)
def test_at_positions(model_string):
    mo = getattr(img, model_string)()
    rng = np.random.default_rng(1)
    x, y, z = rng.uniform(-15., 15., size=(3, 4, 5))
    values = mo.at_positions(x, y, z)
    expected = np.array([[np.asarray(mo.at_position(*p), dtype=float) for p in zip(*row)] for row in zip(x, y, z)])
    expected = np.moveaxis(expected, -1, 0) if expected.ndim == 3 else expected
    assert values.shape == ((3, 4, 5) if model_string in regular_models else (4, 5))
    assert np.array_equal(values, expected, equal_nan=True)


def test_at_positions_broadcasting():
    mo = img.JF12RegularField()
    line = np.linspace(-10., 10., 7)
    values = mo.at_positions(line, 2., 0.1)
    assert values.shape == (3, 7)
    assert np.array_equal(values[:, 3], np.asarray(mo.at_position(0., 2., 0.1)))
    assert mo.at_positions(1., 2., 3.).shape == (3,)
    assert np.array_equal(img.AxiSymmetricSpiral().at_positions(line, 1., 0.)[:, 0], img.AxiSymmetricSpiral().at_position(-10., 1., 0.))
    with pytest.raises(ValueError):
        mo.at_positions(np.zeros(3), np.zeros(4), 0.)


def test_model_variants():
    uf = img.UFMagneticField(model="expX")
    assert uf.model == "expX"
    assert uf.fPoloidalA == uf.all_parameters["expX"]["fPoloidalA"]
    reference = uf.at_position(-8.5, 1., .3)
    uf.set_model("spur")
    uf.set_model("expX")
    assert uf.at_position(-8.5, 1., .3) == reference
    assert img.UFMagneticField().model == "base"
    with pytest.raises(AttributeError):
        uf.model = "base"
    with pytest.raises(ValueError):
        uf.set_model("unknown")

    tf = img.TFMagneticField("Dd1", "C1")
    assert (tf.disk_model, tf.halo_model) == ("Dd1", "C1")
    assert np.isfinite(tf.at_position(5., 3., .5)).all()
    assert img.TFMagneticField().H_disk == 0.055
    with pytest.raises(ValueError):
        tf.set_model("Xd1", "C0")


def test_configuration_members():
    han = img.HanMagneticField()
    position = (-8.5, 1., .1)
    assert np.any(np.asarray(han.at_position(*position)) != 0.)
    han.R_max = 8.
    assert np.all(np.asarray(han.at_position(*position)) == 0.)
    assert len(han.R_s) == 7

    ymw = img.YMW16()
    reference = ymw.at_position(-8.4, .05, .02)
    ymw.localbubble_boundary = .05
    assert ymw.at_position(-8.4, .05, .02) != reference
    for name in ["t0_theta0", "h0", "h1", "h2", "Xgc", "Ygc", "Zgc", "t5_lc", "t6_zyl1", "t6_zyl2", "max_radius"]:
        assert isinstance(getattr(ymw, name), float)

    uf = img.UFMagneticField()
    assert uf.fMaxRadius == 30.
    uf.fMaxRadius = 20.
    assert np.all(np.asarray(uf.at_position(-22., 0., 0.)) == 0.)


def test_point_cloud():
    rng = np.random.default_rng(3)
    positions = rng.uniform(-15., 15., size=(20, 3))
    cloud = img.PointCloud.from_positions(positions)
    assert len(cloud) == cloud.size == 20
    assert np.array_equal(cloud.y, positions[:, 1])
    same = img.PointCloud(positions[:, 0], positions[:, 1], positions[:, 2])
    for model in [img.JF12RegularField(), img.UFMagneticField(), img.YMW16()]:
        on_cloud = model.evaluate(cloud)
        expected = model.at_positions(positions[:, 0], positions[:, 1], positions[:, 2])
        assert on_cloud.shape == expected.shape
        assert np.array_equal(on_cloud, expected)
        assert np.array_equal(model.evaluate(same), on_cloud)

    class Constant(img.RegularVectorField):
        def at_position(self, x, y, z):
            return [x, 2. * y, 3.]

    assert np.array_equal(Constant().evaluate(cloud), np.vstack([positions[:, 0], 2. * positions[:, 1], np.full(20, 3.)]))

    with pytest.raises(ValueError):
        img.PointCloud([1., 2.], [1.], [1., 2.])
    with pytest.raises(ValueError):
        img.PointCloud.from_positions(np.zeros((4, 2)))


def test_jf12_variants():
    jf12 = img.JF12RegularField()
    assert jf12.model == "JF12" and jf12.arm_shift == 1. and jf12.b_arm_6 == -4.2 and jf12.B0_X == 4.6
    reference = jf12.at_position(-8.5, 1., .3)
    jf12.b_arm_2 = 9.
    jf12.set_model("Planck12b")
    assert (jf12.b_arm_2, jf12.b_arm_6, jf12.B0_X, jf12.arm_shift) == (3., -3.5, 1.8, 1.)
    planck_c = img.JF12RegularField("Planck12c")
    assert (planck_c.Bn, planck_c.Bs, planck_c.B0_X, planck_c.b_arm_2, planck_c.b_arm_4, planck_c.b_arm_5, planck_c.b_arm_6,
            planck_c.arm_shift) == (1., -0.8, 3., 2., 2., -3., -3.5, 0.97)
    jf12.set_model("JF12")
    assert jf12.at_position(-8.5, 1., .3) == reference
    with pytest.raises(ValueError):
        jf12.set_model("unknown")
