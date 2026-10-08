# IMAGINE Model Library

A library of Galactic models — regular and random magnetic fields and thermal electron densities — for Galactic inference engines.
The models are written in C++ and can be used from C++ and from Python.
A few (non-essential) models are written in pure Python and are only available from Python.
All implemented models are listed [below](#list-of-models).

> **Work in progress.** Interfaces may still change, and not every model has been validated against its original publication.
> Please report anything that looks wrong.

**Contents**: [Requirements](#requirements) · [Installation](#installation) · [Quick start (Python)](#quick-start-python) · [Quick start (C++)](#quick-start-c) · [Conventions](#conventions) · [Adding a model](#adding-a-model) · [Development](#development) · [List of models](#list-of-models) · [License](#license)


## Requirements

- A C++17 compiler and [CMake](https://cmake.org/) ≥ 3.16
- For the Python package: [Python](https://www.python.org/) ≥ 3.8 and [NumPy](https://numpy.org/) ≥ 1.22

Optional (detected automatically at build time):

- [FFTW3](http://fftw.org/) ≥ 3.3, for the random field models
- [autodiff](https://autodiff.github.io/) (tested with 0.6.12 and 1.1.2) and [Eigen](https://eigen.tuxfamily.org/) ≥ 3.4 (tested with 3.4 and 5.0), for derivatives with respect to model parameters
- [matplotlib](https://matplotlib.org/), for `plot_slice` and the demo notebooks
- [healpy](https://healpy.readthedocs.io/), for the Local Bubble model

[pybind11](https://github.com/pybind/pybind11) is included as a git submodule.


## Installation

### Python

```bash
pip install git+https://github.com/IMAGINE-Consortium/imagine-models.git
```

Append `@branch-name` to the URL to install a specific branch.
The build detects FFTW, autodiff and Eigen if they are installed. To switch features off, set the environment variables `USE_FFTW=OFF` and/or `USE_AUTODIFF=OFF` before running pip.
Check what was built with:

```python
import ImagineModels as img

print(img.has_fftw, img.has_autodiff)
```

To work on the library itself, clone the repository with its submodules and install it in editable mode:

```bash
git clone --recursive https://github.com/IMAGINE-Consortium/imagine-models.git
cd imagine-models
pip install -e ".[test]"
```

Changes to the Python files take effect immediately; after changing C++ code, re-run the `pip install` command.
Alternatively, let the extension rebuild itself on import whenever C++ sources changed (this needs `scikit-build-core`, `cmake` and `ninja` in the environment):

```bash
pip install scikit-build-core cmake ninja
pip install --no-build-isolation -C build-dir=build/python/dev -C editable.rebuild=true -C editable.verbose=false -e ".[test]"
```

### C++

```bash
git clone --recursive https://github.com/IMAGINE-Consortium/imagine-models.git
cd imagine-models
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix /path/to/prefix    # omit --prefix for a system-wide install
```

Optional dependencies are detected automatically; switch them off with `-DUSE_FFTW=OFF` and/or `-DUSE_AUTODIFF=OFF`.
This installs the headers, the library, a CMake package (`ImagineModels::ImagineModels`) and a pkg-config file (`ImagineModels.pc`).


## Quick start (Python)

```python
import numpy as np
import ImagineModels as img

jf12 = img.JF12MagneticField()

# a single position (Galactocentric Cartesian coordinates in kpc), result in µG
jf12.at_position(-8.5, 0.0, 0.1)

# many positions at once, with numpy broadcasting (here: a line along x)
jf12.at_positions(np.linspace(-20.0, 20.0, 401), 0.0, 0.1)  # shape (3, 401)

# on a grid
grid = img.RegularGrid(shape=[200, 200, 40], reference_point=[-20.0, -20.0, -4.0], increment=[0.2, 0.2, 0.2])
b = jf12.evaluate(grid)  # shape (3, 200, 200, 40)
irregular = img.IrregularGrid(x=[-10.0, 2.0], y=[-5.0, 0.0, 4.0], z=[0.0])
img.YMW16().evaluate(irregular)  # scalar field: shape (2, 3, 1)

# on an unstructured set of points (e.g. a catalogue)
cloud = img.PointCloud.from_positions(
    np.array([[-8.5, 0.0, 0.1], [3.0, 4.0, -1.0], [0.0, 12.0, 2.0]])
)  # or PointCloud(x, y, z)
jf12.evaluate(cloud)  # shape (3, 3): (component, point)

# parameters
jf12.parameter_names  # ordered list
jf12.b_arm_1 = 1.2  # single parameter
jf12.parameters = {"Bn": 1.5}  # several at once (partial update)

# published model variants: set_model selects the variant and loads its parameters
uf = img.UF24MagneticField(model="expX")
uf.set_model("spur")

# derivatives w.r.t. the parameters (if built with autodiff)
jf12.active_parameters = ["b_arm_1", "Bn"]
jf12.derivative(-8.5, 1.0, 0.1)  # shape (3, 2)
jf12.derivative(cloud)  # shape (3, 3, 2): (component, point, parameter); also on grids

# random fields (if built with FFTW)
random_field = img.JF12RandomField()
sample = random_field.sample(grid, seed=23)  # shape (3, 200, 200, 40)
random_field.rms(-8.5, 0.0, 0.0)  # analytic local rms amplitude
random_field.rms(grid)  # ... on a grid
img.interpolate(sample, grid, cloud)  # sample interpolated at other positions, shape (3, 3)
```

`demos/python/model_library_demo.ipynb` is a tutorial covering this in more detail, including how to prototype your own model in Python.
`demos/python/model_examples/` has one notebook with plots per model.


## Quick start (C++)

```cpp
#include <iostream>
#include "ImagineModels/ImagineModels.h"

int main() {
    imagine::JF12MagneticField jf12;
    imagine::Vec3<double> b = jf12.at_position(-8.5, 0., 0.1);
    std::cout << b[0] << " " << b[1] << " " << b[2] << std::endl;

    jf12.parameters.b_arm_1 = 1.2;
    imagine::VectorGridData grid_data = jf12.evaluate(imagine::RegularGrid({200, 200, 40}, {-20., -20., -4.}, {.2, .2, .2}));
    // grid_data(component, flat_index), flat_index = (i * ny + j) * nz + k

    imagine::VectorGridData on_points = jf12.evaluate(imagine::PointCloud({-8.5, 3.}, {0., 4.}, {0.1, -1.}));
    // on_points(component, point), shape {2, 1, 1}

#if IMAGINE_HAS_FFTW
    imagine::JF12RandomField random_field;
    imagine::VectorGridData sample = random_field.sample(imagine::RegularGrid({64, 64, 16}, {-10., -10., -2.}, {.3, .3, .25}), 23);
#endif
}
```

With CMake:

```cmake
find_package(ImagineModels REQUIRED)
target_link_libraries(my_target PRIVATE ImagineModels::ImagineModels)
```

or with pkg-config: `c++ -std=c++17 main.cc $(pkg-config --cflags --libs ImagineModels)`.
`IMAGINE_HAS_FFTW` and `IMAGINE_HAS_AUTODIFF` (from `ImagineModels/config.h`) tell which optional parts the installed library contains.
A complete example is `demos/cpp/example.cc`.


## Conventions

- **Units**: positions in kpc, magnetic fields in µG, thermal electron densities in cm⁻³, angles in degrees.
- **Coordinates**: Galactocentric Cartesian (x, y, z), with z perpendicular to the Galactic plane.
- **Grids**: a `RegularGrid` has `shape`, `reference_point` (the first grid point) and `increment`; an `IrregularGrid` has arbitrary x, y and z axes; a `PointCloud` is an unstructured list of N positions. Grid results are indexed `[component, i, j, k]`, point-cloud results `[component, point]` (in C++ the shape is `{N, 1, 1}`). Random fields can be sampled on a `RegularGrid` only (FFT); their `rms` works on all three.
- **Interpolation**: `interpolate(data, grid, points)` (C++: `ImagineModels/Interpolation.h`) evaluates data on a `RegularGrid` (e.g. a random sample) at other positions, linearly or at the nearest grid point; positions outside the grid raise an error unless `nan_outside` is set. Linear interpolation reduces the variance between grid points and does not keep a divergence-free field divergence-free; it is only meaningful if the grid resolves the correlation length of the field.
- **Random fields**: a random field is `rms(x) · G(x)`, where `G` is a zero-mean Gaussian random field with unit variance. The power spectrum only sets the correlation structure and is normalised on the given grid, so a sample realises exactly `rms(x)²` as local variance, independent of the grid. For vector fields `rms` is the total field strength `sqrt(E|B|²)`. Divergence cleaning (`clean_divergence`, on by default) keeps the total power, but with a spatially varying `rms` the local amplitude then follows `rms(x)` only approximately. `GaussianScalarField` is `mu + sigma · G`, `LogNormalScalarField` is `exp(log_mu + log_sigma · G)`.


## Adding a model

A regular model is a C++ class with a single parameter list and a templated field function (`Han.h`/`han.cc` are a compact example):

```cpp
// c_library/headers/ImagineModels/MyModel.h
#define MYMODEL_PARAMETERS(X) \
    X(b0, 2.)                 \
    X(r0, 5.)

IMAGINE_PARAMETERS(MyModelParameters, MYMODEL_PARAMETERS)

class MyMagneticField : public RegularVectorModel<MyMagneticField, MyModelParameters> {
public:
    double r_max = 20.;     // settings that are not fit parameters: plain members

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const MyModelParameters<T> &p) const;
};

// c_library/source/mymodel.cc: implement field<T> (parameters via p.b0, settings via r_max), then
IMAGINE_INSTANTIATE_VECTOR_MODEL(MyMagneticField)
```

Each parameter is declared exactly once; parameter access by name, `derivative` and the Python attributes are generated from that list.
To make the model available, add the source file to `c_library/CMakeLists.txt` and the header to `ImagineModels/RegularModels.h`.
For Python, add a binding file `python_wrapper/regular/mymodel.cc` with `bind_regular_model<MyMagneticField>(m, "MyMagneticField")` (plus `.def_readwrite` for settings), declare it in `python_wrapper/bindings.h`, call it in `python_wrapper/module.cc`, add it to the source list in the top-level `CMakeLists.txt`, and export the class in `ImagineModels/__init__.py`.

Quick prototypes can also be written in pure Python by subclassing `img.RegularVectorField` (see the tutorial notebook).

### Code style

- Formatting is done by tools: `clang-format` (C++, `.clang-format`) and `ruff format` / `ruff check` (Python, `pyproject.toml`). CI checks both with the versions pinned in `.github/workflows/ci.yml`.
- One file stem per model: `ModelName.h`, `model_name.cc`, `python_wrapper/regular/model_name.cc`, `model_name_demo.ipynb`. Headers use `#pragma once`. The Python class name equals the C++ class name.
- Parameter names follow the symbols of the publication or reference code; units other than kpc and µG are given in the parameter list (`X(b_p, -12.) /* deg */`). Angle inputs are in degrees and converted where used (`p.b_p * units::deg`); constants such as `units::pi` come from `ImagineModels/units.h` (no `M_PI`).
- Call math functions on templated values unqualified (`sqrt(x)`, `exp(x)`) so that the autodiff overloads apply.
- Code from other projects is used only under GPL-3.0-compatible licences (GPL, LGPL, BSD, MIT, Apache-2.0, CC-BY-4.0); copyright notices of the original files are kept. Models whose reference code has no licence are written from the publication only.
- Each model header starts with its reference, the external code it is based on (with its licence), and its deviations from the publication (also listed [above](#deviations-from-the-publications)):

  ```cpp
  // Reference: Sun et al. 2008, arXiv:0711.1572
  // Based on: hammurabi v3.01, GPL-3.0
  // Deviations:
  // - halo with the parameters of Sun & Reich 2010
  ```

- Otherwise, comments are at most five words: names of logical blocks, important switches, or equation labels (`// eq. 6`).


## Development

```bash
pytest                                                    # Python tests
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build && ctest --test-dir build   # C++ tests
```

- `tests/test_reference.py` compares every model against stored reference data (`tests/reference_data/`) and checks derivatives against finite differences. If a change alters a model's output on purpose, regenerate the affected cases with `python tests/generate_reference.py --force CASE...`.
- `tests/test_external.py` compares models with values computed by their external reference implementations (`tests/external_data/`, values only): the original YMW16 code, hammurabiX (Jaffe), CRPropa (TF17, Pshirkov, JF12, Archimedes, KST24), the KST24 authors' code, the NE2025/NE2001 Fortran code, the authors' UF23 code and the old hammurabi (Sun disk, Stanev, Fauvet). Each file records the source, tolerance and any deliberate difference.
- `tests/test_random_fields.py` checks the statistics of the random fields (variance, amplitude, divergence, anisotropy).
- The C++ tests (`c_library/test/`, [Catch2](https://github.com/catchorg/Catch2) v3: an installed copy is used if found, otherwise it is downloaded at configure time) check every model for grid consistency, finite values, the parameter registry and derivatives against finite differences, plus the random-field statistics. `build/c_library/imagine_tests --list-tests` lists them; a tag such as `"[random]"` runs a subset.
- CI (`.github/workflows/ci.yml`) checks formatting and lint, runs the C++ and Python tests with all optional dependencies and without them, and builds a wheel from the source distribution.


## List of models

Each model follows its reference publication; the publications are its documentation.
"Based on" names the external code an implementation was ported from or compared with, and its licence.
Differences to the publications are listed under [Deviations from the publications](#deviations-from-the-publications); [CHANGELOG.md](CHANGELOG.md) lists changes between versions.

### Magnetic (vector) fields

| Model | Class | Reference | Based on | Variants (`set_model`) | Notebook |
| ----- | ----- | --------- | -------- | ---------------------- | -------- |
| **Regular models** | | | | | |
| Uniform | `UniformMagneticField` | | | | |
| Helix | `HelixMagneticField` | | | | |
| Archimedean spiral | `ArchimedeanMagneticField` | [Jokipii et al. (1977)](https://ui.adsabs.harvard.edu/abs/1977ApJ...213..861J/abstract) | [CRPropa](https://github.com/CRPropa/CRPropa3) (GPL-3.0) | | [ipynb](demos/python/model_examples/archimedes_demo.ipynb) |
| Axisymmetric spiral (Python only) | `AxiSymmetricSpiral` | | V. Pelgrims | | |
| Local Bubble (Python only)¹ | `LBMagneticField` | [Pelgrims et al. (2020)](https://www.aanda.org/articles/aa/full_html/2020/04/aa37157-19/aa37157-19.html) | V. Pelgrims | | |
| Fauvet | `FauvetMagneticField` | [Fauvet et al. (2012)](https://arxiv.org/abs/1201.5742) | [hammurabi v3.01](https://sourceforge.net/projects/hammurabicode/) (GPL-3.0) | | [ipynb](demos/python/model_examples/fauvet_demo.ipynb) |
| Han | `HanMagneticField` | [Han et al. (2018)](https://iopscience.iop.org/article/10.3847/1538-4365/aa9c45); XH24 disk: [Xu & Han (2024)](https://arxiv.org/abs/2404.02038) | XH24 disk: parameter values from the [authors' code](http://zmtt.bao.ac.cn/GMF/)² | Han2018, XH24 | [ipynb](demos/python/model_examples/han_demo.ipynb) |
| HMR | `HMRMagneticField` | [Harari et al. (1999)](https://arxiv.org/abs/astro-ph/9906309) | [hammurabi v3.01](https://sourceforge.net/projects/hammurabicode/) (GPL-3.0) | | [ipynb](demos/python/model_examples/hmr_demo.ipynb) |
| Jaffe | `JaffeMagneticField` | [Jaffe et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010MNRAS.401.1013J/abstract); Jaffe13: [Jaffe et al. (2013)](https://arxiv.org/abs/1302.0143) | [hammurabiX](https://github.com/hammurabi-dev/hammurabiX) (GPL-3.0); Jaffe13: [hammurabi v3.01](https://sourceforge.net/projects/hammurabicode/) (GPL-3.0) | hammurabiX, Jaffe13 | [ipynb](demos/python/model_examples/jaffe_demo.ipynb) |
| JF12 | `JF12MagneticField` | [Jansson & Farrar (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...757...14J/abstract); Planck variants: [Planck XLII (2016)](https://arxiv.org/abs/1601.00546) | [hammurabiX](https://github.com/hammurabi-dev/hammurabiX) (GPL-3.0), [CRPropa](https://github.com/CRPropa/CRPropa3) (GPL-3.0) | JF12, Planck12b, Planck12c | [ipynb](demos/python/model_examples/jf12_demo.ipynb) |
| KST24 | `KST24MagneticField` | [Korochkin, Semikoz & Tinyakov (2025)](https://arxiv.org/abs/2407.02148) | [CRPropa](https://github.com/CRPropa/CRPropa3) (GPL-3.0); [authors' code](https://doi.org/10.5281/zenodo.14743599) (CC-BY-4.0) | | |
| Pshirkov | `PshirkovMagneticField` | [Pshirkov et al. (2011)](https://iopscience.iop.org/article/10.1088/0004-637X/738/2/192) | [CRPropa](https://github.com/CRPropa/CRPropa3) (GPL-3.0) | ASS, BSS | [ipynb](demos/python/model_examples/pshirkov_demo.ipynb) |
| Stanev | `StanevBSSMagneticField` | [Stanev (1997)](https://arxiv.org/abs/astro-ph/9607086) | [hammurabi v3.01](https://sourceforge.net/projects/hammurabicode/) (GPL-3.0) | | [ipynb](demos/python/model_examples/stanev_demo.ipynb) |
| Sun | `SunMagneticField` | [Sun et al. (2008)](https://www.aanda.org/articles/aa/abs/2008/02/aa8671-07/aa8671-07.html); halo: [Sun & Reich (2010)](https://arxiv.org/abs/1010.4394) | [hammurabi v3.01](https://sourceforge.net/projects/hammurabicode/) (GPL-3.0) | | [ipynb](demos/python/model_examples/sun_demo.ipynb) |
| SVT22 | `SVT22MagneticField` | [Shaw et al. (2022)](https://academic.oup.com/mnras/article/517/2/2534/6731784) | | | [ipynb](demos/python/model_examples/svt22_demo.ipynb) |
| TF17 | `TF17MagneticField` | [Terral & Ferrière (2017)](https://arxiv.org/abs/1611.10222); field forms: [Ferrière & Terral (2014)](https://arxiv.org/abs/1312.1974) (halo model C, disks from A, B, D) | [CRPropa](https://github.com/CRPropa/CRPropa3) (GPL-3.0) | disk Ad1, Bd1, Dd1 × halo C0, C1 | [ipynb](demos/python/model_examples/tf17_demo.ipynb) |
| TT | `TTMagneticField` | [Tinyakov & Tkachev (2002)](https://arxiv.org/abs/astro-ph/0111305) | [hammurabi v3.01](https://sourceforge.net/projects/hammurabicode/) (GPL-3.0) | | [ipynb](demos/python/model_examples/tt_demo.ipynb) |
| UF24 | `UF24MagneticField` | [Unger & Farrar (2024)](https://arxiv.org/abs/2311.12120) | authors' code ([UF23Field v1.1](https://doi.org/10.5281/zenodo.11321212)) (BSD-2-Clause) | base, neCL, expX, spur, cre10, synCG, twistX, nebCor | [ipynb](demos/python/model_examples/uf24_demo.ipynb) |
| WMAP | `WMAPMagneticField` | [Page et al. (2007)](https://iopscience.iop.org/article/10.1086/513699) | [hammurabi v3.01](https://sourceforge.net/projects/hammurabicode/) (GPL-3.0) | | [ipynb](demos/python/model_examples/wmap_demo.ipynb) |
| XH24 halo | `XH24MagneticField` | [Xu & Han (2024)](https://arxiv.org/abs/2404.02038) | compared with the [authors' code](http://zmtt.bao.ac.cn/GMF/)² | | |
| **Random models** | | | | | |
| ES | `ESRandomField` | | [hammurabiX](https://github.com/hammurabi-dev/hammurabiX) (GPL-3.0) | | |
| JF12 | `JF12RandomField` | [Jansson & Farrar (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...761L..11J/abstract); Planck variants: [Planck XLII (2016)](https://arxiv.org/abs/1601.00546) | [hammurabiX](https://github.com/hammurabi-dev/hammurabiX) (GPL-3.0), [CRPropa](https://github.com/CRPropa/CRPropa3) (GPL-3.0) | JF12, Planck12b, Planck12c | [ipynb](demos/python/model_examples/jf12_random_demo.ipynb) |
| UF26 | `UF26RandomField` | [Unger & Farrar (2026)](https://arxiv.org/abs/2608.21293) | | expDisk, ringDisk | |

¹ `from ImagineModels.MagneticFields.LocalBubbleMagneticField import LBMagneticField`; defined on the shell only, requires healpy.

² No licence stated; no code copied.

### Thermal electron (scalar) fields

| Model | Class | Reference | Based on | Variants (`set_model`) | Notebook |
| ----- | ----- | --------- | -------- | ---------------------- | -------- |
| **Regular models** | | | | | |
| Uniform | `UniformDensityField` | | | | |
| NE2025 | `NE2025` | [Ocker & Cordes (2026)](https://arxiv.org/abs/2602.11838); NE2001: [Cordes & Lazio (2002)](https://arxiv.org/abs/astro-ph/0207156) | authors' Fortran code in [mwprop](https://github.com/stella-ocker/mwprop) (GPL-3.0-or-later) | NE2025, NE2001 | |
| Plane-parallel | `PlaneParallelDensity` | [Ocker, Cordes & Chatterjee (2020)](https://arxiv.org/abs/2004.11921); Gaensler08: [Gaensler et al. (2008)](https://arxiv.org/abs/0808.2550) | | Ocker20, Gaensler08 | |
| YMW16 | `YMW16` | [Yao et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...835...29Y/abstract) | original C code v1.3.1 (via [pygedm](https://github.com/FRBs/pygedm)) (GPL-3.0-or-later) | | [ipynb](demos/python/model_examples/ymw16_demo.ipynb) |
| YT20 halo | `YT20` | [Yamasaki & Totani (2020)](https://arxiv.org/abs/1909.00849) | | | |
| **Random models** | | | | | |
| Gaussian | `GaussianScalarField` | | | | [ipynb](demos/python/model_examples/gaussian_scalar_demo.ipynb) |
| Log-normal | `LogNormalScalarField` | | | | |

### Deviations from the publications

Models not listed here have no known deviations.

- **Archimedean spiral**: no fitted model; dimensionless parameters as in CRPropa (`R_0` in kpc, `Omega / v_w` in 1/kpc, `B_0` the radial field at `R_0`).
- **ES random field**: no publication; rms profile as in hammurabiX, `b0 * sqrt(exp(-(r - r_obs)/r0) exp(-(|z| - |z_obs|)/z0))`.
- **Han, variant XH24**: disk as in the authors' code, not in the papers: `R_s(6)` = 8.16 kpc, an extra zone 10.5–15 kpc with `B_s7` = 4.5 µG, disk to 20 kpc.
- **Jaffe**: 3D form and default parameters from the hammurabiX template, not from a publication (the 2010 model is 2D, with R1 = 3 kpc and an arm cutoff at 15 kpc).
- **Jaffe, variant Jaffe13**: Table A1 azimuths read as clockwise, so the arms are at 350°, 260°, 170°, 80° (amplitudes 3, 0.5, −4, 1.2 µG); this matches the NE2001 arms and puts the −4 µG arm at the Sagittarius-Carina arm. Arm height h_c = 2 kpc from Table A1 (Planck XLII quotes 0.5 kpc as the original value). Field directions, arm distances and the reversal inside a negative ring as in hammurabi v3.01; no reversal above the disk (introduced only for Jaffe13b).
- **JF12 (regular)**: `b8` from flux conservation (2.755 µG; the paper rounds to 2.7). Molecular ring field `b_ring · 5 kpc / r` as in CRPropa and hammurabiX (the paper gives no radial dependence). Planck variants: striation factor β not modelled, `b8` from flux conservation (CRPropa keeps 2.7).
- **JF12 (random)**: Planck variants without the striation factor β.
- **KST24**: from the authors' code, not in the paper: Sagittarius-Carina arm widening by 3° along the arm, radial (3–17 kpc) and vertical arm cut-offs, arm widths capped at 1.2 kpc, spiral scale a = 3 kpc; Sagittarius-Carina `rdisk_sagcar` = 0.79 kpc as in the code (Table 2: 0.8); outer Perseus field −3.5 µG as in Table 2 and CRPropa (the authors' Zenodo code uses −2.5 µG).
- **NE2025 / NE2001**: electron density only (the fluctuation parameters F and scattering are not included); double instead of single precision; position parameters (Galactic Centre, local ISM) in the NE2001 frame: x towards l = 90°, Sun at (0, 8.5, 0) kpc.
- **Plane-parallel, Ocker20**: smooth plane-parallel component only; the paper's clumps and voids are specific to single lines of sight.
- **Stanev**: eq. 4 used with exp(−|z|/z0) (sign missing in the paper); field cut at cylindrical r = 20 kpc (the paper: 20 kpc in all directions).
- **Sun**: halo with the parameters of Sun & Reich (2010): `bH_B0` = 2 µG, `bH_z1a`/`bH_z1b` = 0.2/4 kpc.
- **SVT22**: `B_val` = 3.72 µG; the published best fit is 3.96 µG (3 µG in arXiv v1).
- **TF17**: lower limits of Table 2 used as values, as in CRPropa.
- **UF26**: the paper constrains only the rms; the power spectrum is the library default.
- **WMAP**: ψ0 = 27° and sin ψ on r̂ as corrected in Jansson et al. (2009, Sec. 5.3.6); the paper gives no amplitude, `b_b0` = 6 µG is a default; `anti` (field reversed for z > 0) from hammurabi, not in the paper.
- **XH24 halo**: field set to zero beyond r = 20 kpc, as in the authors' code.
- **YT20**: physical constants at full precision (the authors' script and pygedm use three digits), Υ = 2.61; density set to zero beyond r_vir, where the paper ends the DM integration.
- **YMW16**: Galactic part only (no Magellanic Clouds or IGM), R ≤ 30 kpc as in the original code; exact π instead of the original's `RAD = 57.295779` (relative differences ≤ 4e-6), Gum Nebula limit θ → 0 instead of 0/0, no gap for azimuths in [6.28, 2π).


## License

GPL-3.0-or-later (see [COPYING](COPYING)).
Models adapted from other codes are listed with their origin and its licence in the [model tables](#list-of-models); all are compatible with GPL-3.0, and copyright notices of the original files are kept in the corresponding source files.
