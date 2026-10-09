"""Local Bubble magnetic field models (Python only).

Alves et al. 2018 (thin spheroidal shell), Pelgrims et al. 2020 (thin shell from 3D dust maps),
Pelgrims et al. 2025 (thick shell) and O'Neill et al. 2024 (data). Positions are Galactocentric
in kpc with the Sun at ``sun_position``; angles are in degrees.
"""

import numpy as np

from ImagineModels import RegularVectorField

SUN_POSITION = (-8.5, 0.0, 0.0)  # kpc


def _healpy():
    try:
        import healpy
    except ImportError as error:
        raise ImportError("This Local Bubble model requires healpy (pip install healpy).") from error
    return healpy


def _direction(lon, lat):
    lon, lat = np.radians(lon), np.radians(lat)
    return np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def _unit(v):
    return v / np.linalg.norm(v, axis=0)


def _dot(a, b):
    return np.sum(a * b, axis=0)


def _thin_shell_direction(normal, e_r, b0):
    # A18 eq. 6
    b0 = b0[:, None]
    return _unit(b0 * _dot(normal, e_r) - e_r * _dot(normal, b0))


def _tangents(e):
    polar = np.abs(e[2]) > 0.9
    helper = np.array([polar, np.zeros_like(polar), ~polar], dtype=float)
    t1 = _unit(np.cross(e, helper, axis=0))
    return t1, np.cross(e, t1, axis=0)


class _ThinShellField(RegularVectorField):
    """Direction of the field on a thin shell, evaluated along lines of sight from the Sun."""

    def __init__(self):
        super().__init__()
        self.sun_position = SUN_POSITION

    def at_LonLat(self, lon, lat):
        """Unit field direction (3, n) on the shell toward Galactic (lon, lat) [deg]."""
        raise NotImplementedError

    def position_at_LonLat(self, lon, lat):
        """Galactocentric position (3, n) [kpc] of the shell toward (lon, lat) [deg]."""
        raise NotImplementedError

    def at_position(self, x, y, z):
        raise ValueError(
            "This model is defined on the shell only; use at_LonLat() for the field and position_at_LonLat() "
            "for the shell positions."
        )


class Alves18MagneticField(_ThinShellField):
    """Alves et al. 2018 (A&A 611, L5): field direction in a thin spheroidal shell, Table 1."""

    def __init__(self):
        super().__init__()
        self.l0 = 71.0  # deg
        self.b0 = -16.0  # deg
        self.a_ell = 0.1  # kpc, equatorial semi-axis
        self.c_ratio = 2.7  # c_ell / a_ell
        self.psi_ell = 216.0  # deg, precession
        self.theta_ell = 30.0  # deg, nutation
        self.delta = (0.17, 0.56, 0.0)  # units of a_ell

    def _geometry(self):
        center = np.asarray(self.delta, dtype=float) * self.a_ell
        psi, theta = np.radians(self.psi_ell), np.radians(self.theta_ell)
        axis = np.array([np.sin(psi) * np.sin(theta), -np.cos(psi) * np.sin(theta), np.cos(theta)])  # ZXZ
        c_ell = self.c_ratio * self.a_ell

        def scale(v):
            return v / self.a_ell + (1.0 / c_ell - 1.0 / self.a_ell) * _dot(v, axis[:, None]) * axis[:, None]

        return center, scale

    def _shell(self, lon, lat):
        u = _direction(np.atleast_1d(lon), np.atleast_1d(lat))
        center, scale = self._geometry()
        su, sc = scale(u), scale(np.broadcast_to(center[:, None], u.shape))
        a, b, c = _dot(su, su), _dot(su, sc), _dot(sc, sc) - 1.0
        distance = (b + np.sqrt(b * b - a * c)) / a
        point = distance * u
        normal = _unit(scale(scale(point - center[:, None])))
        return point, normal, center

    def at_LonLat(self, lon, lat):
        point, normal, center = self._shell(lon, lat)
        return _thin_shell_direction(normal, _unit(point - center[:, None]), _direction(self.l0, self.b0))

    def position_at_LonLat(self, lon, lat):
        point, _, _ = self._shell(lon, lat)
        return point + np.asarray(self.sun_position)[:, None]


class Pelgrims20MagneticField(_ThinShellField):
    """Pelgrims et al. 2020 (A&A 636, A17): field direction in the thin shell of the Local Bubble.

    ``shell_map``: path to ``L19_map-inner_final.fits`` (https://doi.org/10.7910/DVN/RHPVNC) or a HEALPix
    (RING) map of heliocentric distances to the shell in kpc.
    """

    available_models = ("lmax2", "lmax4", "lmax6", "lmax8", "lmax10")
    _table1 = {  # Table 1
        "lmax2": (32.0, 15.4, -170.6, 14.9, 71.6),
        "lmax4": (-16.9, -184.5, -195.5, 15.8, 73.1),
        "lmax6": (57.6, 79.2, -86.3, 16.8, 73.2),
        "lmax8": (-9.0, -96.7, -150.2, 14.3, 72.9),
        "lmax10": (51.2, 121.3, -107.6, 13.0, 72.6),
    }

    def __init__(self, shell_map, model="lmax6"):
        super().__init__()
        self._shell_map = shell_map
        self.set_model(model)

    def set_model(self, model):
        """Load the Table 1 parameters and, for a file, the shell surface of the given l_max."""
        if model not in self.available_models:
            raise ValueError(f"Unknown Pelgrims20 model '{model}'.")
        self.model = model
        dx, dy, dz, self.b0, self.l0 = self._table1[model]
        self.dx, self.dy, self.dz = dx / 1000.0, dy / 1000.0, dz / 1000.0  # kpc
        if isinstance(self._shell_map, str):
            field = 1 + self.available_models.index(model)
            distance = _healpy().read_map(self._shell_map, field=field) / 1000.0
        else:
            distance = np.asarray(self._shell_map, dtype=float)
        self.update_shell(distance)

    def update_shell(self, distance):
        """Set the shell from a HEALPix map of heliocentric distances [kpc]."""
        hp = _healpy()
        self.nside = hp.get_nside(distance)
        npix = hp.nside2npix(self.nside)
        self.surface = distance
        self.shell_points = np.asarray(hp.pix2vec(self.nside, np.arange(npix))) * distance
        # normal from neighbouring pixels
        neighbours = hp.get_all_neighbours(self.nside, np.arange(npix))[np.arange(0, 8, 2), :]
        p = self.shell_points[:, neighbours]
        normal = _unit(np.cross((p[:, 2] - p[:, 0]).T, (p[:, 1] - p[:, 3]).T).T)
        self.normal = normal * np.sign(_dot(_unit(self.shell_points), normal))

    def _pixels(self, lon, lat):
        return _healpy().ang2pix(self.nside, np.atleast_1d(lon), np.atleast_1d(lat), lonlat=True)

    def at_LonLat(self, lon, lat):
        pix = self._pixels(lon, lat)
        point = self.shell_points[:, pix]
        center = np.array([[self.dx], [self.dy], [self.dz]])
        return _thin_shell_direction(self.normal[:, pix], _unit(point - center), _direction(self.l0, self.b0))

    def position_at_LonLat(self, lon, lat):
        return self.shell_points[:, self._pixels(lon, lat)] + np.asarray(self.sun_position)[:, None]


class ONeill24MagneticField(_ThinShellField):
    """O'Neill et al. 2024 (ApJ): 3D field orientation on the Local Bubble surface (data model).

    ``table``: path to ``ONeill2024_LocalBubble_3DBfield.fits`` (https://doi.org/10.7910/DVN/A8HWUF). The
    orientation is a pseudovector (B and -B are equivalent); NaN where the model is undefined.
    """

    nside = 256

    def __init__(self, table):
        super().__init__()
        try:
            from astropy.io import fits
        except ImportError as error:
            raise ImportError("ONeill24MagneticField requires astropy (pip install astropy).") from error
        with fits.open(table) as hdul:
            data = hdul[1].data
            self.shell_points = np.array([data["x"], data["y"], data["z"]]) / 1000.0  # kpc, heliocentric
            self.orientation = np.array([data["Bx"], data["By"], data["Bz"]])
            self.delta_theta = np.array(data["delta_thetaB"])  # deg

    def _pixels(self, lon, lat):
        return _healpy().ang2pix(self.nside, np.atleast_1d(lon), np.atleast_1d(lat), lonlat=True)

    def at_LonLat(self, lon, lat):
        return self.orientation[:, self._pixels(lon, lat)]

    def position_at_LonLat(self, lon, lat):
        return self.shell_points[:, self._pixels(lon, lat)] + np.asarray(self.sun_position)[:, None]

    def uncertainty_at_LonLat(self, lon, lat):
        """Orientation uncertainty delta_thetaB [deg]."""
        return self.delta_theta[self._pixels(lon, lat)]


class _Sphere:
    def __init__(self, center, radius):
        self.center, self.radius = np.asarray(center, dtype=float), radius

    def from_center(self, c, e):
        c = np.asarray(c, dtype=float)
        d = (c - self.center)[:, None]
        de = _dot(d, e)
        r = -de + np.sqrt(de * de - _dot(d, d) + self.radius**2)
        n = (c[:, None] + r * e - self.center[:, None]) / self.radius
        return r, r * (e - n / _dot(n, e))


class _RadialMap:
    """Closed surface as a HEALPix map of distances from an origin (heliocentric, kpc)."""

    def __init__(self, origin, distance):
        self.origin, self.distance = np.asarray(origin, dtype=float), np.asarray(distance, dtype=float)
        hp = _healpy()
        self.nside = hp.get_nside(self.distance)
        self.step = hp.nside2resol(self.nside)

    def radius(self, u):
        hp = _healpy()
        theta = np.arccos(np.clip(u[2], -1.0, 1.0))
        return hp.get_interp_val(self.distance, theta, np.arctan2(u[1], u[0]))

    def ray(self, c, e, iterations=60):
        """Distance from point c along unit vectors e to the surface (bisection; c inside)."""
        c = np.asarray(c, dtype=float)[:, None]
        lo = np.zeros(e.shape[1])
        hi = np.full(e.shape[1], 2.0 * self.distance.max() + 2.0 * np.linalg.norm(c[:, 0] - self.origin))
        for _ in range(iterations):
            mid = 0.5 * (lo + hi)
            q = c + mid * e - self.origin[:, None]
            outside = np.linalg.norm(q, axis=0) > self.radius(_unit(q))
            hi, lo = np.where(outside, mid, hi), np.where(outside, lo, mid)
        return 0.5 * (lo + hi)

    def seen_from(self, c):
        hp = _healpy()
        e = np.asarray(hp.pix2vec(self.nside, np.arange(self.distance.size)))
        return _RadialMap(c, self.ray(c, e))

    def from_center(self, c, e):
        table = self if np.allclose(self.origin, c) else self.seen_from(c)
        r = table.radius(e)
        t1, t2 = _tangents(e)
        h = table.step
        gradient = 0.0
        for t in (t1, t2):
            plus = table.radius(np.cos(h) * e + np.sin(h) * t)
            minus = table.radius(np.cos(h) * e - np.sin(h) * t)
            gradient = gradient + (plus - minus) / (2.0 * h) * t
        return r, gradient


class _CachedSurface:
    """Surface tabulated about the explosion centre for fast evaluation."""

    def __init__(self, surface, center):
        self.center = np.asarray(center, dtype=float)
        if isinstance(surface, _Sphere):
            self._sphere, self._table = surface, None
        else:
            same = np.allclose(surface.origin, center)
            self._sphere, self._table = None, surface if same else surface.seen_from(center)

    def __call__(self, e):
        if self._sphere is not None:
            return self._sphere.from_center(self.center, e)
        return self._table.from_center(self.center, e)


class Pelgrims25MagneticField(RegularVectorField):
    """Pelgrims, Unger & Maris 2025 (A&A 695, A148): magnetic field in the thick shell of the Local Bubble.

    Scenarios (Table 1): SCO, SCA (spherical shells), DCO, DCA (inner surface from ``shell_map``, constant
    thickness), DDO, DDA (inner and outer surfaces from ``shell_map`` and ``outer_map``). ``shell_map``: path
    to ``L19_map-inner_final.fits`` (https://doi.org/10.7910/DVN/RHPVNC, column l_max = 6) or a HEALPix map
    of heliocentric distances [kpc]; ``outer_map``: HEALPix map [kpc] of the outer surface (not public).
    """

    available_models = ("SCO", "SCA", "DCO", "DCA", "DDO", "DDA")

    def __init__(self, model="SCO", shell_map=None, outer_map=None):
        super().__init__()
        self.sun_position = SUN_POSITION
        self.shell_only = True  # zero outside the shell
        self._shell_map, self._outer_map = shell_map, outer_map
        self.set_model(model)

    def set_model(self, model):
        """Select a scenario of Table 1 and load its parameters."""
        if model not in self.available_models:
            raise ValueError(f"Unknown Pelgrims25 model '{model}'.")
        self.model = model
        # Table 1, heliocentric [kpc]
        self.center_sphere = (-0.0248, -0.0326, -0.0233)
        self.radius_sphere = 0.2167
        self.thickness = 0.035
        self.center_p20 = (0.023, -0.034, -0.122)
        self.l0, self.b0 = 73.0, 17.0  # deg
        self.B0 = 3.0  # muG, Table 2
        self.update_shell()

    def update_shell(self):
        """Rebuild the shell surfaces after changing parameters."""
        center = self.center_sphere if self.model[2] == "O" else self.center_p20
        self.explosion_center = np.asarray(center, dtype=float)
        o = np.asarray(self.center_sphere, dtype=float)
        if self.model.startswith("SC"):
            inner = _Sphere(o, self.radius_sphere)
            outer = _Sphere(o, self.radius_sphere + self.thickness)
        else:
            if self._shell_map is None:
                raise ValueError(f"Scenario {self.model} needs shell_map (inner surface of Pelgrims et al. 2020).")
            if isinstance(self._shell_map, str):
                distance = _healpy().read_map(self._shell_map, field=3) / 1000.0  # l_max = 6
            else:
                distance = np.asarray(self._shell_map, dtype=float)
            inner = _RadialMap((0.0, 0.0, 0.0), distance)
            if self.model.startswith("DC"):
                inner_o = inner.seen_from(o)
                outer = _RadialMap(o, inner_o.distance + self.thickness)
            else:
                if self._outer_map is None:
                    raise ValueError(f"Scenario {self.model} needs outer_map (outer surface, Appendix B).")
                outer = _RadialMap((0.0, 0.0, 0.0), np.asarray(self._outer_map, dtype=float))
        self._inner = _CachedSurface(inner, self.explosion_center)
        self._outer = _CachedSurface(outer, self.explosion_center)

    def _field(self, x, y, z):
        p = np.array([x, y, z], dtype=float) - np.asarray(self.sun_position)[:, None]
        q = p - self.explosion_center[:, None]
        r = np.linalg.norm(q, axis=0)
        e = np.where(r > 0.0, q / np.where(r > 0.0, r, 1.0), np.array([[0.0], [0.0], [1.0]]))
        r_min, d_min = self._inner(e)
        r_max, d_max = self._outer(e)
        b0 = self.B0 * _direction(self.l0, self.b0)[:, None]
        b0_r = _dot(b0, e)
        b0_t = b0 - b0_r * e
        width = r_max - r_min
        in_shell = (r >= r_min) & (r <= r_max)
        r_safe = np.where(in_shell, r, 1.0)
        ratio = r_max * (r - r_min) / (width * r_safe)  # r0 / r
        grad = (r_safe * (r_min * d_max - r_max * d_min) + r_max**2 * d_min - r_min**2 * d_max) / width**2
        b = ratio**2 * b0_r * e + ratio * _dot(grad, b0_t) / r_safe * e + ratio * (r_max / width) * b0_t  # eq. 13
        outside = 0.0 if self.shell_only else b0
        return np.where(in_shell, b, np.where(r > r_max, outside, 0.0))

    def at_position(self, x, y, z):
        return self._field([x], [y], [z])[:, 0].tolist()

    def at_positions(self, x, y, z):
        """Vectorised evaluation, broadcasting x, y, z; returns shape (3, ...)."""
        x, y, z = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (x, y, z)))
        return self._field(x.ravel(), y.ravel(), z.ravel()).reshape((3,) + x.shape)
