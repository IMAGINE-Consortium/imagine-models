import ImagineModels as img

import pytest
import numpy as np


if not img.__has_random_fields__:
    pytest.skip("ImagineModels was built without FFTW", allow_module_level=True)

vector_models = ['JF12RandomField', 'ESRandomField', 'UF26RandomField']
scalar_models = ['GaussianScalarField', 'LogNormalScalarField']

grid = img.RegularGrid(shape=[16, 16, 16], reference_point=[-1., -1., -1.], increment=[.125, .125, .125])
stat_grid = img.RegularGrid(shape=[48, 32, 24], reference_point=[-8., -4., -2.], increment=[.25, .25, .25])
seeds = range(12)


def _sample(model_string, seed):
    return getattr(img, model_string)().sample(grid, seed)


def _within(values, expected, n_sigma=5.):
    values = np.asarray(values)
    error = values.std(ddof=1) / np.sqrt(len(values))
    assert abs(values.mean() - expected) < n_sigma * max(error, 1e-12), (values.mean(), expected, error)


class ConstantRandomField(img.RandomVectorField):
    def __init__(self, rms=2., slope=None):
        super().__init__()
        self.value = rms
        self.slope = slope
        self.apply_spectrum = slope is not None

    def rms(self, x, y, z):
        return self.value

    def spectrum(self, k):
        return k ** -self.slope


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
@pytest.mark.parametrize('apply_spectrum', [False, True])
def test_unit_variance(model_string, apply_spectrum):
    model = getattr(img, model_string)()
    model.apply_spectrum = apply_spectrum
    variances = []
    for seed in seeds:
        g = model.random_numbers(stat_grid, seed)
        assert abs(g.mean()) < 1e-10 if g.ndim == 3 else np.allclose(g.mean(axis=(1, 2, 3)), 0., atol=1e-10)
        variances.append((g ** 2).sum(axis=0).mean() if g.ndim == 4 else (g ** 2).mean())
    _within(variances, 1.)


def test_gaussian_mean_and_rms():
    model = img.GaussianScalarField()
    model.mu, model.sigma, model.apply_spectrum = 2.5, 0.3, False
    samples = np.array([model.sample(stat_grid, seed) for seed in seeds])
    _within(samples.mean(axis=(1, 2, 3)), 2.5)
    _within(samples.var(axis=(1, 2, 3)), 0.3 ** 2)
    assert model.mean(0., 0., 0.) == 2.5
    assert model.rms(0., 0., 0.) == 0.3


def test_lognormal_mean_and_rms():
    model = img.LogNormalScalarField()
    model.log_mu, model.log_sigma, model.apply_spectrum = 0.2, 0.5, False
    samples = np.array([model.sample(stat_grid, seed) for seed in seeds])
    assert model.mean(0., 0., 0.) == pytest.approx(np.exp(0.2 + 0.125))
    assert model.rms(0., 0., 0.) == pytest.approx(np.sqrt(np.expm1(0.25)) * np.exp(0.2 + 0.125))
    _within(samples.mean(axis=(1, 2, 3)), float(model.mean(0., 0., 0.)))
    _within(samples.var(axis=(1, 2, 3)), float(model.variance(0., 0., 0.)))


@pytest.mark.parametrize('slope', [None, 2.])
def test_vector_amplitude_matches_rms(slope):
    model = ConstantRandomField(rms=2., slope=slope)
    model.clean_divergence = False
    energies = [(model.sample(stat_grid, seed) ** 2).sum(axis=0).mean() for seed in seeds]
    _within(energies, 4.)


@pytest.mark.parametrize('model_string', vector_models)
def test_model_amplitude_follows_rms_profile(model_string):
    model = getattr(img, model_string)()
    model.clean_divergence = False
    model.apply_spectrum = False
    rms = model.rms(stat_grid)
    selected = rms > 1e-3
    ratios = [((model.sample(stat_grid, seed) ** 2).sum(axis=0)[selected] / rms[selected] ** 2).mean() for seed in seeds]
    _within(ratios, 1.)


@pytest.mark.parametrize('model_string', vector_models + scalar_models)
def test_rms_interfaces(model_string):
    model = getattr(img, model_string)()
    irregular = img.IrregularGrid(np.asarray([-8.5, 3.]), np.asarray([0., 4.]), np.asarray([0.1]))
    on_grid = model.rms(irregular)
    assert on_grid.shape == (2, 2, 1)
    assert on_grid[1, 1, 0] == model.rms(3., 4., .1)
    line = model.rms(np.array([-8.5, 3.]), 4., .1)
    assert line.shape == (2,)
    assert line[1] == on_grid[1, 1, 0]
    assert np.allclose(model.variance(np.array([-8.5, 3.]), 4., .1), line ** 2)
    cloud = img.PointCloud([-8.5, 3.], [4., 4.], [.1, .1])
    assert np.array_equal(model.rms(cloud), line)


def test_python_subclass():
    model = ConstantRandomField(rms=1.5)
    assert model.rms(1., 2., 3.) == 1.5
    assert model.sample(grid, 1).shape == (3, *grid.shape)


def test_sample_requires_regular_grid():
    irregular = img.IrregularGrid(np.asarray([0., 1.]), np.asarray([2.]), np.asarray([-1., 0., 1.]))
    with pytest.raises(TypeError):
        img.GaussianScalarField().sample(irregular, 3)


def _relative_divergence(field, grid):
    b = np.fft.fftn(field, axes=(1, 2, 3))
    k = np.meshgrid(*[np.fft.fftfreq(n, d) for n, d in zip(grid.shape, grid.increment)], indexing='ij')
    divergence = k[0] * b[0] + k[1] * b[1] + k[2] * b[2]
    scale = np.sqrt(k[0] ** 2 + k[1] ** 2 + k[2] ** 2) * np.sqrt((np.abs(b) ** 2).sum(axis=0))
    return np.abs(divergence).sum() / scale.sum()


@pytest.mark.parametrize('model', [ConstantRandomField(rms=2.), ConstantRandomField(rms=2., slope=2.), img.JF12RandomField(), img.ESRandomField()])
def test_divergence_cleaning(model):
    model.clean_divergence = True
    assert _relative_divergence(model.sample(stat_grid, 4), stat_grid) < 1e-12
    model.clean_divergence = False
    assert _relative_divergence(model.sample(stat_grid, 4), stat_grid) > 0.3


@pytest.mark.parametrize('slope', [None, 2.])
def test_divergence_cleaning_preserves_amplitude(slope):
    model = ConstantRandomField(rms=2., slope=slope)
    energies = [(model.sample(stat_grid, seed) ** 2).sum(axis=0).mean() for seed in seeds]
    _within(energies, 4.)


@pytest.mark.parametrize('model_string', vector_models)
def test_divergence_cleaning_preserves_total_power(model_string):
    model = getattr(img, model_string)()
    rms2 = (model.rms(stat_grid) ** 2).sum()
    ratios = [(model.sample(stat_grid, seed) ** 2).sum() / rms2 for seed in seeds]
    _within(ratios, 1.)


class VerticalAnisotropy(ConstantRandomField):
    def anisotropy_direction(self, x, y, z):
        return [0., 0., 3.]


@pytest.mark.parametrize('rho', [1., 2., 0.5])
def test_anisotropy(rho):
    model = VerticalAnisotropy(rms=2.)
    model.clean_divergence = False
    model.anisotropy_rho = rho
    samples = [model.sample(stat_grid, seed) for seed in seeds]
    _within([(b ** 2).sum(axis=0).mean() for b in samples], 4.)
    _within([(b[2] ** 2).mean() / (b ** 2).sum(axis=0).mean() for b in samples], rho ** 4 / (rho ** 4 + 2.))
    model.apply_anisotropy = False
    _within([(b[2] ** 2).mean() / (b ** 2).sum(axis=0).mean() for b in [model.sample(stat_grid, seed) for seed in seeds]], 1. / 3.)


def test_jf12_anisotropy_follows_regular_field():
    model = img.JF12RandomField()
    regular = img.JF12RegularField()
    assert np.allclose(model.anisotropy_direction(-8.5, 1., .2), regular.at_position(-8.5, 1., .2))


def _vertical_integral(model, r, z_max=30., n=60001):
    z = np.linspace(0., z_max, n)
    b2 = model.rms(np.full(n, -r), 0., z) ** 2
    return 0.5 * np.sum(b2[1:] + b2[:-1]) * (z[1] - z[0])


def test_uf26_against_paper():
    disk = img.UF26RandomField()
    assert disk.model == "expDisk"
    assert disk.b_disk == 4.4 and disk.z_disk == 1.0 and disk.l_r == 15. and disk.r_c == 1.
    sun = disk.r_sun
    assert disk.rms(-sun, 0., disk.z_disk) / disk.rms(-sun, 0., 0.) == pytest.approx(np.exp(-1.))
    assert disk.rms(-sun, 0., 0.) == pytest.approx(4.4 / (1. + np.exp((sun - 18.) / 2.)))
    assert _vertical_integral(disk, sun) == pytest.approx(11., abs=1.)

    ring = img.UF26RandomField("ringDisk")
    assert (ring.b_disk, ring.z_disk, ring.l_r, ring.b_ring, ring.r_ring, ring.w_ring, ring.z_ring) == (4.1, 0.9, 100., 4.1, 4.5, 2.0, 1.9)
    assert _vertical_integral(ring, sun) == pytest.approx(10., abs=1.)
    ring_only = img.UF26RandomField("ringDisk")
    ring_only.b_disk = 0.
    assert _vertical_integral(ring_only, ring.r_ring) == pytest.approx(30., abs=9.)

    disk.set_model("ringDisk")
    assert np.array_equal(disk.rms(np.array([-8., 3.]), 1., .2), ring.rms(np.array([-8., 3.]), 1., .2))
    with pytest.raises(ValueError):
        img.UF26RandomField("unknown")
