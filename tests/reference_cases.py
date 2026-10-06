import re
from pathlib import Path

import numpy as np

import ImagineModels as img


REPO = Path(__file__).resolve().parents[1]
DATA_DIR = Path(__file__).resolve().parent / "reference_data"

RTOL_REGULAR, ATOL_REGULAR = 1e-12, 1e-14
RTOL_RANDOM, ATOL_RANDOM = 1e-9, 1e-12

KNOWN_BAD_JACOBIAN = {"HanMagneticField": "B3", "StanevBSSMagneticField": "B4", "SunMagneticField": "B5"}
NO_GRID_CASES = {}
BROKEN_DERIVATIVE_CASES = {"TFMagneticField__Dd1_C0": "B15", "TFMagneticField__Dd1_C1": "B15",
                           "UFMagneticField__expX": "B16", "SVT22__default": "B17"}


def _positions():
    rng = np.random.default_rng(20261006)
    special = np.array([[0., 0., 0.], [-8.5, 0., 0.], [8.5, 0., 0.], [0., 8.5, 0.], [-8.5, 0., 1.],
                        [-8.5, 0., -1.], [3., 4., 0.2], [-1., -1., -1.], [12., -9., 2.], [0., 0., 5.]])
    box = rng.uniform([-20., -20., -5.], [20., 20., 5.], size=(54, 3))
    return np.vstack([special, box])


POSITIONS = _positions()
JACOBIAN_POSITIONS = POSITIONS[[1, 4, 6, 8, 11, 17, 23, 42]]

REGULAR_GRID = dict(shape=[5, 4, 3], reference_point=[-15., -12., -3.], increment=[6., 6., 2.5])
IRREGULAR_GRID = dict(grid_x=np.array([-17., -8.5, -2., 0.5, 6., 14.]),
                      grid_y=np.array([-11., 0., 3.5, 9.]),
                      grid_z=np.array([-2., 0., 0.4, 3.]))

RANDOM_GRIDS = {"even": dict(shape=[8, 8, 8], reference_point=[-4., -4., -4.], increment=[1., 1., 1.]),
                "odd": dict(shape=[7, 6, 5], reference_point=[-3., -2., -1.], increment=[.5, .7, .9])}
RANDOM_SEEDS = [3, 7]


def _regular_cases():
    cases = {}
    for name in ["ArchimedeanMagneticField", "FauvetMagneticField", "HMRMagneticField", "HanMagneticField",
                 "HelixMagneticField", "JF12RegularField", "JaffeMagneticField", "PshirkovMagneticField",
                 "SVT22", "StanevBSSMagneticField", "SunMagneticField", "TFMagneticField", "TTMagneticField",
                 "UFMagneticField", "UniformDensityField", "UniformMagneticField", "WMAPMagneticField", "YMW16"]:
        cases[f"{name}__default"] = (name, {})
    cases["UniformMagneticField__set"] = ("UniformMagneticField", {"bx": -3.2, "by": 1.5, "bz": 0.25})
    cases["UniformDensityField__set"] = ("UniformDensityField", {"n0": 0.03})
    cases["JF12RegularField__no_halo"] = ("JF12RegularField", {"do_halo": False})
    cases["JF12RegularField__no_X"] = ("JF12RegularField", {"do_X": False})
    cases["JaffeMagneticField__ring_no_bar"] = ("JaffeMagneticField", {"ring": True, "bar": False})
    cases["JaffeMagneticField__bss"] = ("JaffeMagneticField", {"bss": True})
    cases["JaffeMagneticField__quadruple"] = ("JaffeMagneticField", {"quadruple": True})
    cases["PshirkovMagneticField__ass"] = ("PshirkovMagneticField", {"useASS": True, "useBSS": False})
    cases["PshirkovMagneticField__no_halo"] = ("PshirkovMagneticField", {"useHalo": False})
    cases["WMAPMagneticField__anti"] = ("WMAPMagneticField", {"b_anti": True})
    for disk in ["Ad1", "Bd1", "Dd1"]:
        for halo in ["C0", "C1"]:
            cases[f"TFMagneticField__{disk}_{halo}"] = ("TFMagneticField", {"activeDiskModel": disk, "activeHaloModel": halo})
    for variant in ["base", "neCL", "expX", "spur", "cre10", "synCG", "twistX", "nebCor"]:
        cases[f"UFMagneticField__{variant}"] = ("UFMagneticField", {"activeModel": variant, "set_parameters": variant})
    cases["AxiSymmetricSpiral__default"] = ("AxiSymmetricSpiral", {})
    return cases


def _random_cases():
    cases = {}
    for name in ["JF12RandomField", "ESRandomField", "GaussianScalarField", "LogNormalScalarField"]:
        cases[f"{name}__default"] = (name, {})
    cases["JF12RandomField__no_cleaning"] = ("JF12RandomField", {"clean_divergence": False})
    cases["JF12RandomField__no_spectrum"] = ("JF12RandomField", {"apply_spectrum": False})
    return cases


REGULAR_CASES = _regular_cases()
RANDOM_CASES = _random_cases() if img.__has_random_fields__ else {}


def make_model(name, settings):
    model = getattr(img, name)()
    for key, value in settings.items():
        if key == "set_parameters":
            model.set_parameters(value)
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
    jac = np.atleast_2d(np.asarray(model.derivative(*position), dtype=float))
    if hasattr(model, "parameter_names"):
        names = list(model.parameter_names)
        jac = jac[:, [names.index(c) for c in columns]]
    return jac


def sample(model, grid, seed):
    return np.asarray(model.sample(img.RegularGrid(**grid), seed), dtype=float)


def parameter_defaults(model):
    out = {}
    for key, prop in type(model).__dict__.items():
        if isinstance(prop, property) and prop.fset is not None:
            value = getattr(model, key)
            if isinstance(value, (bool, int, float, str)) or (isinstance(value, list) and all(isinstance(v, (int, float)) for v in value)):
                out[key] = value
    return out


_SOURCE_FILES = {"ArchimedeanMagneticField": "archimedes.cc", "FauvetMagneticField": "fauvet.cc",
                 "HMRMagneticField": "hararimollerachroulet.cc", "HanMagneticField": "han.cc",
                 "HelixMagneticField": "helix.cc", "JF12RegularField": "regularjf12.cc",
                 "JaffeMagneticField": "jaffe.cc", "PshirkovMagneticField": "pshirkov.cc", "SVT22": "svt22.cc",
                 "StanevBSSMagneticField": "stanevbss.cc", "SunMagneticField": "sun.cc",
                 "TFMagneticField": "tf17.cc", "TTMagneticField": "tinyakovtkachev.cc",
                 "UFMagneticField": "ungerfarrar.cc", "WMAPMagneticField": "wmap.cc", "YMW16": "ymw16.cc"}


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
