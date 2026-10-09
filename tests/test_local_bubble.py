import numpy as np
import pytest

from ImagineModels.MagneticFields.LocalBubbleMagneticField import (
    Alves18MagneticField,
    ONeill24MagneticField,
    Pelgrims20MagneticField,
    Pelgrims25MagneticField,
)

SUN = np.array([[-8.5], [0.0], [0.0]])


def _lonlat(v):
    v = v / np.linalg.norm(v)
    return np.degrees(np.arctan2(v[1], v[0])) % 360.0, np.degrees(np.arcsin(v[2]))


def _sphere_distance(center, radius, u):
    d = -np.asarray(center)[:, None]
    du = np.sum(d * u, axis=0)
    return -du + np.sqrt(du**2 - np.sum(d * d, axis=0) + radius**2)


def _shell_points(model, n=300, seed=2, margin=0.002):
    rng = np.random.default_rng(seed)
    e = rng.normal(size=(3, n))
    e /= np.linalg.norm(e, axis=0)
    r = rng.uniform(model.radius_sphere + margin, model.radius_sphere + model.thickness - margin, n)
    return np.asarray(model.center_sphere)[:, None] + e * r + SUN


def test_alves18_against_paper():
    model = Alves18MagneticField()
    lon = np.linspace(0.0, 350.0, 36).repeat(12)
    lat = np.tile(np.linspace(-85.0, 85.0, 12), 36)
    b = model.at_LonLat(lon, lat)
    assert np.allclose(np.linalg.norm(b, axis=0), 1.0)
    point, normal, center = model._shell(lon, lat)
    assert np.abs(np.sum(b * normal, axis=0)).max() < 1e-12
    # field in the plane of B0 and e_r
    b0 = np.array([np.cos(np.radians(-16)) * np.cos(np.radians(71)), np.cos(np.radians(-16)) * np.sin(np.radians(71))])
    plane = np.cross(np.append(b0, np.sin(np.radians(-16))), (point - center[:, None]).T)
    assert np.abs(np.sum(plane.T * b, axis=0)).max() < 1e-12
    # Table 1 fit: mean directions over the polar caps (|b| > 60 deg), paper (70, 43) and (74, -14)
    grid_lon, grid_lat = np.meshgrid(np.arange(0.5, 360.0, 1.0), np.arange(-89.5, 90.0, 1.0))
    weight = np.cos(np.radians(grid_lat))
    for cap, expected in ((grid_lat >= 60, (70.0, 43.0)), (grid_lat <= -60, (74.0, -14.0))):
        mean = (model.at_LonLat(grid_lon[cap], grid_lat[cap]) * weight[cap]).sum(axis=1)
        lon_b, lat_b = _lonlat(mean)
        assert abs(lon_b - expected[0]) < 11 and abs(lat_b - expected[1]) < 8
    position = model.position_at_LonLat(lon, lat)
    assert np.allclose(position - SUN, point)
    with pytest.raises(ValueError):
        model.at_position(-8.5, 0.0, 0.0)


def test_pelgrims25_spherical_bubble():
    model = Pelgrims25MagneticField("SCO")
    assert (model.B0, model.l0, model.b0, model.thickness) == (3.0, 73.0, 17.0, 0.035)
    model.center_sphere = (0.0, 0.0, 0.0)
    model.update_shell()
    r_min, r_max = model.radius_sphere, model.radius_sphere + model.thickness
    rng = np.random.default_rng(1)
    e = rng.normal(size=(3, 200))
    e /= np.linalg.norm(e, axis=0)
    r = rng.uniform(r_min, r_max, 200)
    b = model.at_positions(*(e * r + SUN))
    b0 = np.array([np.cos(np.radians(17)) * np.cos(np.radians(73)), np.cos(np.radians(17)) * np.sin(np.radians(73))])
    cos_a = np.sum(np.append(b0, np.sin(np.radians(17)))[:, None] * e, axis=0)
    s = (r - r_min) / r
    expected = 3.0 * s * (r_max / (r_max - r_min)) ** 2 * np.sqrt(1 - cos_a**2 + s**2 * cos_a**2)  # eq. 21
    assert np.allclose(np.linalg.norm(b, axis=0), expected, rtol=1e-10)
    assert np.allclose(model.at_position(-8.5, 0.0, 0.0), 0.0)
    assert np.allclose(model.at_positions(-8.5 + 0.3, 0.0, 0.0), 0.0)
    model.shell_only = False
    assert np.allclose(model.at_positions(-8.5 + 0.3, 0.0, 0.0), 3.0 * np.append(b0, np.sin(np.radians(17))))


def test_pelgrims25_divergence_free():
    model = Pelgrims25MagneticField("SCA")
    p = _shell_points(model)
    h = 1e-5
    div = sum(
        (
            model.at_positions(*(p + h * np.eye(3)[i][:, None]))[i]
            - model.at_positions(*(p - h * np.eye(3)[i][:, None]))[i]
        )
        / (2 * h)
        for i in range(3)
    )
    b = np.linalg.norm(model.at_positions(*p), axis=0)
    assert np.max(np.abs(div) * model.thickness / b) < 1e-5
    assert np.allclose(model.at_positions(*p)[:, 0], model.at_position(*p[:, 0]))
    with pytest.raises(ValueError):
        Pelgrims25MagneticField("DCO")
    with pytest.raises(ValueError):
        Pelgrims25MagneticField("unknown")


def test_pelgrims25_map_surfaces():
    hp = pytest.importorskip("healpy")
    sco = Pelgrims25MagneticField("SCA")
    nside = 64
    u = np.asarray(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))))
    inner = _sphere_distance(sco.center_sphere, sco.radius_sphere, u)
    outer = _sphere_distance(sco.center_sphere, sco.radius_sphere + sco.thickness, u)
    p = _shell_points(sco, n=200, margin=0.005)
    expected = sco.at_positions(*p)
    for model in (Pelgrims25MagneticField("DCA", shell_map=inner), Pelgrims25MagneticField("DDA", inner, outer)):
        diff = np.linalg.norm(model.at_positions(*p) - expected, axis=0) / np.linalg.norm(expected, axis=0)
        assert np.median(diff) < 1e-3 and diff.max() < 2e-2


def test_pelgrims20_thin_shell():
    hp = pytest.importorskip("healpy")
    nside = 32
    u = np.asarray(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))))
    model = Pelgrims20MagneticField(np.full(u.shape[1], 0.2), model="lmax2")
    assert (model.dx, model.l0, model.b0) == (0.032, 71.6, 14.9)
    model.set_model("lmax6")
    assert np.allclose((model.dx, model.dy, model.dz, model.l0, model.b0), (0.0576, 0.0792, -0.0863, 73.2, 16.8))
    lon, lat = hp.pix2ang(nside, np.arange(u.shape[1]), lonlat=True)
    b = model.at_LonLat(lon, lat)
    # sphere around the Sun: normal is radial
    assert np.abs(np.sum(b * u, axis=0)).max() < 0.05
    e_r = u * 0.2 - np.array([[model.dx], [model.dy], [model.dz]])
    b0 = np.array(
        [np.cos(np.radians(16.8)) * np.cos(np.radians(73.2)), np.cos(np.radians(16.8)) * np.sin(np.radians(73.2))]
    )
    plane = np.cross(np.append(b0, np.sin(np.radians(16.8))), e_r.T).T
    assert np.abs(np.sum(plane * b, axis=0)).max() < 0.05
    assert np.allclose(model.position_at_LonLat(lon, lat) - SUN, 0.2 * u)
    with pytest.raises(ValueError):
        model.set_model("lmax3")


def test_oneill24_table(tmp_path):
    hp = pytest.importorskip("healpy")
    fits = pytest.importorskip("astropy.io.fits")
    npix = hp.nside2npix(256)
    u = np.asarray(hp.pix2vec(256, np.arange(npix)))
    columns = [fits.Column(name=n, format="D", array=a) for n, a in zip(("x", "y", "z"), 150.0 * u)]
    columns += [fits.Column(name=n, format="D", array=a) for n, a in zip(("Bx", "By", "Bz"), np.roll(u, 1, axis=0))]
    columns.append(fits.Column(name="delta_thetaB", format="D", array=np.full(npix, 10.0)))
    path = tmp_path / "oneill.fits"
    fits.BinTableHDU.from_columns(columns).writeto(path)
    model = ONeill24MagneticField(str(path))
    lon, lat = np.array([10.0, 200.0]), np.array([30.0, -45.0])
    pix = hp.ang2pix(256, lon, lat, lonlat=True)
    assert np.allclose(model.at_LonLat(lon, lat), np.roll(u, 1, axis=0)[:, pix])
    assert np.allclose(model.position_at_LonLat(lon, lat), 0.15 * u[:, pix] + SUN)
    assert np.allclose(model.uncertainty_at_LonLat(lon, lat), 10.0)
