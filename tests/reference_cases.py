import re
from pathlib import Path

import numpy as np

import ImagineModels as img

REPO = Path(__file__).resolve().parents[1]
DATA_DIR = Path(__file__).resolve().parent / "reference_data"

RTOL_REGULAR, ATOL_REGULAR = 1e-12, 1e-14
RTOL_RANDOM, ATOL_RANDOM = 1e-9, 1e-12

KNOWN_BAD_JACOBIAN = {}
NO_GRID_CASES = {}
BROKEN_DERIVATIVE_CASES = {}


def _positions():
    rng = np.random.default_rng(20261006)
    special = np.array(
        [
            [0.0, 0.0, 0.0],
            [-8.5, 0.0, 0.0],
            [8.5, 0.0, 0.0],
            [0.0, 8.5, 0.0],
            [-8.5, 0.0, 1.0],
            [-8.5, 0.0, -1.0],
            [3.0, 4.0, 0.2],
            [-1.0, -1.0, -1.0],
            [12.0, -9.0, 2.0],
            [0.0, 0.0, 5.0],
        ]
    )
    box = rng.uniform([-20.0, -20.0, -5.0], [20.0, 20.0, 5.0], size=(54, 3))
    return np.vstack([special, box])


POSITIONS = _positions()
JACOBIAN_POSITIONS = POSITIONS[[1, 4, 6, 8, 11, 17, 23, 42]]

REGULAR_GRID = dict(shape=[5, 4, 3], reference_point=[-15.0, -12.0, -3.0], increment=[6.0, 6.0, 2.5])
IRREGULAR_GRID = dict(
    grid_x=np.array([-17.0, -8.5, -2.0, 0.5, 6.0, 14.0]),
    grid_y=np.array([-11.0, 0.0, 3.5, 9.0]),
    grid_z=np.array([-2.0, 0.0, 0.4, 3.0]),
)

RANDOM_GRIDS = {
    "even": dict(shape=[8, 8, 8], reference_point=[-4.0, -4.0, -4.0], increment=[1.0, 1.0, 1.0]),
    "odd": dict(shape=[7, 6, 5], reference_point=[-3.0, -2.0, -1.0], increment=[0.5, 0.7, 0.9]),
}
RANDOM_SEEDS = [3, 7]


def _regular_cases():
    cases = {}
    for name in [
        "ArchimedeanMagneticField",
        "FauvetMagneticField",
        "HMRMagneticField",
        "HanMagneticField",
        "HelixMagneticField",
        "JF12MagneticField",
        "JaffeMagneticField",
        "PshirkovMagneticField",
        "SVT22MagneticField",
        "StanevBSSMagneticField",
        "SunMagneticField",
        "TFMagneticField",
        "TTMagneticField",
        "UFMagneticField",
        "UniformDensityField",
        "UniformMagneticField",
        "WMAPMagneticField",
        "XH24MagneticField",
        "YMW16",
    ]:
        cases[f"{name}__default"] = (name, {})
    cases["UniformMagneticField__set"] = ("UniformMagneticField", {"bx": -3.2, "by": 1.5, "bz": 0.25})
    cases["UniformDensityField__set"] = ("UniformDensityField", {"n0": 0.03})
    cases["JF12MagneticField__no_halo"] = ("JF12MagneticField", {"do_halo": False})
    cases["JF12MagneticField__no_X"] = ("JF12MagneticField", {"do_X": False})
    cases["HanMagneticField__XH24"] = ("HanMagneticField", {"set_model": ["XH24"]})
    cases["JF12MagneticField__Planck12b"] = ("JF12MagneticField", {"set_model": ["Planck12b"]})
    cases["JF12MagneticField__Planck12c"] = ("JF12MagneticField", {"set_model": ["Planck12c"]})
    cases["SVT22MagneticField__no_halo"] = ("SVT22MagneticField", {"do_halo": False})
    cases["JaffeMagneticField__ring_no_bar"] = ("JaffeMagneticField", {"ring": True, "bar": False})
    cases["JaffeMagneticField__bss"] = ("JaffeMagneticField", {"bss": True})
    cases["JaffeMagneticField__quadruple"] = ("JaffeMagneticField", {"quadruple": True})
    cases["PshirkovMagneticField__ass"] = ("PshirkovMagneticField", {"set_model": ["ASS"]})
    cases["PshirkovMagneticField__no_halo"] = ("PshirkovMagneticField", {"useHalo": False})
    cases["WMAPMagneticField__anti"] = ("WMAPMagneticField", {"b_anti": True})
    for disk in ["Ad1", "Bd1", "Dd1"]:
        for halo in ["C0", "C1"]:
            cases[f"TFMagneticField__{disk}_{halo}"] = ("TFMagneticField", {"set_model": [disk, halo]})
    for variant in ["base", "neCL", "expX", "spur", "cre10", "synCG", "twistX", "nebCor"]:
        cases[f"UFMagneticField__{variant}"] = ("UFMagneticField", {"set_model": [variant]})
    cases["AxiSymmetricSpiral__default"] = ("AxiSymmetricSpiral", {})
    return cases


def _random_cases():
    cases = {}
    for name in ["JF12RandomField", "ESRandomField", "UF26RandomField", "GaussianScalarField", "LogNormalScalarField"]:
        cases[f"{name}__default"] = (name, {})
    cases["UF26RandomField__ringDisk"] = ("UF26RandomField", {"set_model": ["ringDisk"]})
    cases["JF12RandomField__Planck12c"] = ("JF12RandomField", {"set_model": ["Planck12c"]})
    cases["JF12RandomField__no_cleaning"] = ("JF12RandomField", {"clean_divergence": False})
    cases["JF12RandomField__no_spectrum"] = ("JF12RandomField", {"apply_spectrum": False})
    return cases


REGULAR_CASES = _regular_cases()
RANDOM_CASES = _random_cases() if img.__has_random_fields__ else {}


RENAMED = {"JF12RegularField": "JF12MagneticField", "SVT22": "SVT22MagneticField"}


def model_class(name):
    return getattr(img, RENAMED.get(name, name), None)


def make_model(name, settings):
    model = model_class(name)()
    for key, value in settings.items():
        if key == "set_model":
            model.set_model(*value)
        else:
            setattr(model, key, value)
    return model


def at_positions(model, positions):
    return np.array([np.atleast_1d(np.asarray(model.at_position(*p), dtype=float)) for p in positions])


def on_regular_grid(model):
    return np.asarray(model.evaluate(img.RegularGrid(**REGULAR_GRID)), dtype=float)


def on_irregular_grid(model):
    g = IRREGULAR_GRID
    return np.asarray(model.evaluate(img.IrregularGrid(g["grid_x"], g["grid_y"], g["grid_z"])), dtype=float)


def has_jacobian(model):
    return img.__has_autodiff__ and hasattr(model, "derivative")


def jacobian(model, position, columns):
    model.active_parameters = list(columns)
    return np.atleast_2d(np.asarray(model.derivative(*position), dtype=float))


def sample(model, grid, seed):
    return np.asarray(model.sample(img.RegularGrid(**grid), seed), dtype=float)


def parameter_defaults(model):
    out = {}
    for key, prop in type(model).__dict__.items():
        if isinstance(prop, property) and prop.fset is not None:
            value = getattr(model, key)
            if isinstance(value, (bool, int, float, str)) or (
                isinstance(value, list) and all(isinstance(v, (int, float)) for v in value)
            ):
                out[key] = value
    return out


_SOURCE_FILES = {
    "ArchimedeanMagneticField": "archimedes.cc",
    "FauvetMagneticField": "fauvet.cc",
    "HMRMagneticField": "hmr.cc",
    "HanMagneticField": "han.cc",
    "HelixMagneticField": "helix.cc",
    "JF12MagneticField": "jf12.cc",
    "JaffeMagneticField": "jaffe.cc",
    "PshirkovMagneticField": "pshirkov.cc",
    "SVT22MagneticField": "svt22.cc",
    "StanevBSSMagneticField": "stanev.cc",
    "SunMagneticField": "sun.cc",
    "TFMagneticField": "tf17.cc",
    "TTMagneticField": "tt.cc",
    "UFMagneticField": "uf24.cc",
    "WMAPMagneticField": "wmap.cc",
    "YMW16": "ymw16.cc",
}


def jacobian_column_labels(name, model):
    if hasattr(model, "parameter_names"):
        return list(model.parameter_names)
    if name == "UniformMagneticField":
        return ["bx", "by", "bz"]
    if name == "UniformDensityField":
        return ["n0"]
    src = (REPO / "c_library" / "source" / _SOURCE_FILES[name]).read_text()
    args = re.search(r"ad::wrt\(([^)]*)\)", src).group(1)
    return [a.strip().split(".")[-1] for a in args.split(",")]
