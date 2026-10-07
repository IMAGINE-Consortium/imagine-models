# IMAGINE Model Library

A library of Galactic models — regular and random magnetic fields and thermal electron densities — for Galactic inference engines.
The models are written in C++ and can be used from C++ and from Python.
A few (non-essential) models are written in pure Python and are only available from Python.
All implemented models are listed [below](#list-of-models).

> **Work in progress.** Interfaces may still change, and not every model has been validated against its original publication.
> Please report anything that looks wrong.

**Contents**: [Requirements](#requirements) · [Installation](#installation) · [Quick start (Python)](#quick-start-python) · [Quick start (C++)](#quick-start-c) · [Conventions](#conventions) · [Adding a model](#adding-a-model) · [Development](#development) · [List of models](#list-of-models)


## Requirements

- A C++17 compiler and [CMake](https://cmake.org/) ≥ 3.16
- For the Python package: [Python](https://www.python.org/) ≥ 3.8 and [NumPy](https://numpy.org/) ≥ 1.22

Optional (detected automatically at build time):

- [FFTW3](http://fftw.org/) ≥ 3.3, for the random field models
- [autodiff](https://autodiff.github.io/) (tested with 0.6.12 and 1.1.2) and [Eigen3](https://eigen.tuxfamily.org/) ≥ 3.4, for derivatives with respect to model parameters
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

jf12 = img.JF12RegularField()

# a single position (Galactocentric Cartesian coordinates in kpc), result in µG
jf12.at_position(-8.5, 0., 0.1)

# many positions at once, with numpy broadcasting (here: a line along x)
jf12.at_positions(np.linspace(-20., 20., 401), 0., 0.1)        # shape (3, 401)

# on a grid
grid = img.RegularGrid(shape=[200, 200, 40], reference_point=[-20., -20., -4.], increment=[.2, .2, .2])
b = jf12.evaluate(grid)                                          # shape (3, 200, 200, 40)
irregular = img.IrregularGrid(x=[-10., 2.], y=[-5., 0., 4.], z=[0.])
img.YMW16().evaluate(irregular)                                  # scalar field: shape (2, 3, 1)

# on an unstructured set of points (e.g. a catalogue)
cloud = img.PointCloud.from_positions(np.array([[-8.5, 0., 0.1], [3., 4., -1.], [0., 12., 2.]]))   # or PointCloud(x, y, z)
jf12.evaluate(cloud)                                             # shape (3, 3): (component, point)

# parameters
jf12.parameter_names                  # ordered list
jf12.b_arm_1 = 1.2                    # single parameter
jf12.parameters = {"Bn": 1.5}         # several at once (partial update)

# published model variants: set_model selects the variant and loads its parameters
uf = img.UFMagneticField(model="expX")
uf.set_model("spur")

# derivatives w.r.t. the parameters (if built with autodiff)
jf12.active_parameters = ["b_arm_1", "Bn"]
jf12.derivative(-8.5, 1., 0.1)        # shape (3, 2)

# random fields (if built with FFTW)
random_field = img.JF12RandomField()
sample = random_field.sample(grid, seed=23)       # shape (3, 200, 200, 40)
random_field.rms(-8.5, 0., 0.)                    # analytic local rms amplitude
random_field.rms(grid)                            # ... on a grid
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

- **Units**: positions in kpc, magnetic fields in µG, thermal electron densities in cm⁻³.
- **Coordinates**: Galactocentric Cartesian (x, y, z), with z perpendicular to the Galactic plane.
- **Grids**: a `RegularGrid` has `shape`, `reference_point` (the first grid point) and `increment`; an `IrregularGrid` has arbitrary x, y and z axes; a `PointCloud` is an unstructured list of N positions. Grid results are indexed `[component, i, j, k]`, point-cloud results `[component, point]` (in C++ the shape is `{N, 1, 1}`). Random fields can be sampled on a `RegularGrid` only (FFT); their `rms` works on all three.
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


## Development

```bash
pytest                                                    # Python tests
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build && ctest --test-dir build   # C++ tests
```

- `tests/test_reference.py` compares every model against stored reference data (`tests/reference_data/`) and checks derivatives against finite differences. If a change alters a model's output on purpose, regenerate the affected cases with `python tests/generate_reference.py --force CASE...`.
- `tests/test_external.py` compares models with values computed by their external reference implementations (`tests/external_data/`, values only): the original YMW16 code, hammurabiX (Jaffe), CRPropa (TF17, Pshirkov, JF12, Archimedes), the authors' UF23 code and the old hammurabi (Sun disk, Stanev, Fauvet). Each file records the source, tolerance and any deliberate difference.
- `tests/test_random_fields.py` checks the statistics of the random fields (variance, amplitude, divergence, anisotropy).
- The C++ tests (`c_library/test/`, [Catch2](https://github.com/catchorg/Catch2) v3: an installed copy is used if found, otherwise it is downloaded at configure time) check every model for grid consistency, finite values, the parameter registry and derivatives against finite differences, plus the random-field statistics. `build/c_library/imagine_tests --list-tests` lists them; a tag such as `"[random]"` runs a subset.
- CI (`.github/workflows/ci.yml`) runs the C++ and Python tests with all optional dependencies and without them, and builds a wheel from the source distribution.


## List of models

### Magnetic (vector) fields

| Model | Python class | C++ | Reference | Notes | Original implementation | Notebook |
| ----- | ------------ | --- | --------- | ----- | ----------------------- | -------- |
| **Regular models** | | | | | | |
| Uniform | `UniformMagneticField` | &#x2714; | | used for unit tests | | |
| Helix | `HelixMagneticField` | &#x2714; | | | | |
| Axisymmetric spiral | `AxiSymmetricSpiral` | &#x2718; | Pelgrims, V. | pure Python | Pelgrims, V. | |
| Archimedean spiral | `ArchimedeanMagneticField` | &#x2714; | | simple demonstrative ASS model; dimensionless as in CRPropa: with `R_0` in kpc, `Omega / v_w` is in 1/kpc and `B_0` is the radial field strength at `R_0` | [CRPropa](https://github.com/CRPropa/CRPropa3) | [ipynb](demos/python/model_examples/archimedes_demo.ipynb) |
| Local Bubble | `LBMagneticField` | &#x2718; | [Pelgrims et al. (2020)](https://www.aanda.org/articles/aa/full_html/2020/04/aa37157-19/aa37157-19.html) | pure Python (`from ImagineModels.MagneticFields.LocalBubbleMagneticField import LBMagneticField`), only defined on the shell, requires healpy | Pelgrims, V. | |
| Jaffe | `JaffeMagneticField` | &#x2714; | [Jaffe et al. (2010)](https://ui.adsabs.harvard.edu/abs/2010MNRAS.401.1013J/abstract) | based on ASS-A spiral with modifications, parameter values taken from hammurabi, not from any publication | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX) | [ipynb](demos/python/model_examples/jaffe_demo.ipynb) |
| Sun2008 | `SunMagneticField` | &#x2714; | [Sun et al. (2008)](https://www.aanda.org/articles/aa/abs/2008/02/aa8671-07/aa8671-07.html) | ASS+Ring as disk field, toroidal asymmetric halo with updated halo parameters from [Sun et al. (2010)](https://iopscience.iop.org/article/10.1088/1674-4527/10/12/009), central part of disk field is constant in z-direction (unphysical) | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/) | [ipynb](demos/python/model_examples/sun_demo.ipynb) |
| Han2018 | `HanMagneticField` | &#x2714; | [Han et al. (2018)](https://iopscience.iop.org/article/10.3847/1538-4365/aa9c45) | BSS-S disk field | | [ipynb](demos/python/model_examples/han_demo.ipynb) |
| Pshirkov | `PshirkovMagneticField` | &#x2714; | [Pshirkov et al. (2011)](https://iopscience.iop.org/article/10.1088/0004-637X/738/2/192) | ASS-S or BSS-S (`set_model("ASS")`/`set_model("BSS")`, default BSS; loads the published pitch angle and southern halo amplitude), disk and halo can be switched off (`useDisk`, `useHalo`); the halo is asymmetric w.r.t. the plane | [CRPropa](https://github.com/CRPropa/CRPropa3) | [ipynb](demos/python/model_examples/pshirkov_demo.ipynb) |
| HMR | `HMRMagneticField` | &#x2714; | [Harari et al. (1999)](https://arxiv.org/abs/astro-ph/9906309) | BSS-S model | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/) and [Kachelrieß (2007)](https://arxiv.org/pdf/astro-ph/0510444.pdf) | [ipynb](demos/python/model_examples/hmr_demo.ipynb) |
| TT | `TTMagneticField` | &#x2714; | [Tinyakov and Tkachev (2002)](https://arxiv.org/abs/astro-ph/0111305) | BSS-A model (eq. 5 in ref.) | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/) and [Kachelrieß (2007)](https://arxiv.org/pdf/astro-ph/0510444.pdf) | [ipynb](demos/python/model_examples/tt_demo.ipynb) |
| TF17 | `TFMagneticField` | &#x2714; | [Terral and Ferrière (2017)](https://arxiv.org/abs/1611.10222) | disk models Ad1/Bd1/Dd1 and halo models C0/C1 (`set_model(disk, halo)`). Only the halo was fitted to data, which leads to very strong field strengths and unexpected features in the disk fields; the halo fields can diverge at large r/z. Better viewed as a mathematical exercise than a model for, e.g., cosmic-ray propagation. | [CRPropa](https://github.com/CRPropa/CRPropa3) | [ipynb](demos/python/model_examples/tf17_demo.ipynb) |
| Fauvet | `FauvetMagneticField` | &#x2714; | [Fauvet et al. (2012)](https://arxiv.org/abs/1201.5742) | modified logarithmic spiral (MLS) with z-component and exponential radial profile `B0 exp(-(r - R_sun)/R_B)` (Sec. 2.1, no halo), pitch angle −30° as used for the simulations in the paper | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/) | [ipynb](demos/python/model_examples/fauvet_demo.ipynb) |
| Stanev | `StanevBSSMagneticField` | &#x2714; | [Stanev (1997)](https://arxiv.org/abs/astro-ph/9607086) | BSS-S model; the change in the halo field at \|z\| = 0.5 kpc was not in the hammurabi implementation | [Hammurabi (old)](https://sourceforge.net/projects/hammurabicode/) | [ipynb](demos/python/model_examples/stanev_demo.ipynb) |
| WMAP | `WMAPMagneticField` | &#x2714; | [Page et al. (2007)](https://iopscience.iop.org/article/10.1086/513699) | logarithmic spiral with constant amplitude and z-component; parameters from the original publication, not from the update in [Ruiz-Granados et al. (2010)](https://www.aanda.org/articles/aa/full_html/2010/14/aa12733-09/aa12733-09.html) | | [ipynb](demos/python/model_examples/wmap_demo.ipynb) |
| Jansson Farrar | `JF12RegularField` | &#x2714; | [Jansson & Farrar (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...757...14J/abstract) | regular JF12 field (disk + symmetric toroidal halo + X-field); `set_model("JF12" | "Planck12b" | "Planck12c")` selects the original or the Planck 2016 updates ([Planck XLII](https://arxiv.org/abs/1601.00546), "Jansson12b/c"; the striated-field factor β is not modelled) | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX) | [ipynb](demos/python/model_examples/jf12_regular_demo.ipynb) |
| Unger Farrar | `UFMagneticField` | &#x2714; | [Unger & Farrar (2024)](https://arxiv.org/abs/2311.12120) | variants base, neCL, expX, spur, cre10, synCG, twistX, nebCor (`set_model(name)`) | Unger & Farrar (BSD-2), see also [CRPropa](https://github.com/CRPropa/CRPropa3) | [ipynb](demos/python/model_examples/uf24_regular_demo.ipynb) |
| SVT22 | `SVT22` | &#x2714; | [Shaw et al. (2022)](https://academic.oup.com/mnras/article/517/2/2534/6731784) | model for the Galactic halo bubble | | [ipynb](demos/python/model_examples/svt22_demo.ipynb) |
| XH24 halo | `XH24MagneticField` | &#x2714; | [Xu & Han (2024)](https://arxiv.org/abs/2404.02038) | toroidal halo field ("huge magnetic toroids", eq. 2, Table 2), antisymmetric w.r.t. the plane; halo only, to be combined with a disk model such as Han2018 | [authors' code](http://zmtt.bao.ac.cn/GMF/) | |
| **Random models** | | | | | | |
| Jansson Farrar | `JF12RandomField` | &#x2714; | [Jansson & Farrar (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...761L..11J/abstract) | rms profile from JF12; anisotropy along the regular JF12 field (`anisotropy_rho`); `set_model("JF12" | "Planck12b" | "Planck12c")` as for the regular field (Planck: B_iso = 7.8 µG with relative arm, interior and halo strengths) | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX) | [ipynb](demos/python/model_examples/jf12_random_demo.ipynb) |
| Ensslin Steininger | `ESRandomField` | &#x2714; | | rms `b0 * sqrt(exp(-(r - r_obs)/r0) exp(-(|z| - |z_obs|)/z0))` (energy density scaled exponentially, `b0` = rms at the observer), as in hammurabiX | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX) | |
| Unger Farrar 2026 | `UF26RandomField` | &#x2714; | [Unger & Farrar (2026)](https://arxiv.org/abs/2608.21293) | isotropic random field, rms profile only (Sec. 8, Table 2): `set_model("expDisk")` (default; sech disk) or `set_model("ringDisk")` (disk + inner annulus, added in quadrature). The paper constrains only the rms; the power spectrum is the library default (`spectral_offset`, `spectral_slope`) | | |

### Thermal electron (scalar) fields

| Model | Python class | C++ | Reference | Notes | Original implementation | Notebook |
| ----- | ------------ | --- | --------- | ----- | ----------------------- | -------- |
| **Regular models** | | | | | | |
| Uniform | `UniformDensityField` | &#x2714; | | used for unit tests | | |
| YMW16 | `YMW16` | &#x2714; | [Yao et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017ApJ...835...29Y/abstract) | | [Hammurabi X](https://github.com/hammurabi-dev/hammurabiX) | [ipynb](demos/python/model_examples/ymw16_demo.ipynb) |
| **Random models** | | | | | | |
| Gaussian | `GaussianScalarField` | &#x2714; | | `mu + sigma · G` | | [ipynb](demos/python/model_examples/gaussian_scalar_demo.ipynb) |
| Log-normal | `LogNormalScalarField` | &#x2714; | | `exp(log_mu + log_sigma · G)` | | |
