# Changelog

## 0.2.0 (unreleased)

Results of many models changed because implementation errors were fixed. Results computed with earlier versions should be checked.

### Changed results

- **YMW16**: now follows the original C code v1.3.1. The outer-disk cutoff had its arguments swapped (densities beyond 15 kpc were too high), the Local arm never contributed, arm widths lacked the factor cos(pitch), the Local Bubble used a wrong scaling and centre, θ_sg was 78.8° instead of 75.8°. Also the warp now acts on the disk components only, and the model is limited to R ≤ 30 kpc, as in the original.
- **Jaffe**: `bar_phi0` is now interpreted in degrees, the field direction inside the bar is fixed, and the peak term vanishes for `r_peak = 0`.
- **TF17**: C1 halos use `phi_star_halo` (was `phi_star_disk`); the Dd1 disk no longer returns NaN; the default disk height is 0.055 kpc (was 0.0055).
- **UF24**: the radial profiles of the `base` and `expX` poloidal halos were swapped; `expX` uses its published `fPoloidalA`; the field extends to 30 kpc (was 20).
- **Pshirkov**: the ASS variant uses its published pitch angle and southern halo amplitude (it used the BSS values).
- **Sun**: ASS+RING disk with R0 = 10 kpc, Rc = 5 kpc; halo antisymmetric at all radii, counter-clockwise in the north, with the Sun & Reich (2010) parameters.
- **WMAP**: ψ0 = 27°; the `anti` switch now has an effect.
- **HMR, TT**: the field follows its own spiral pattern (azimuth convention fixed); TT also has the published phase and halo parity.
- **Fauvet**: now the model of Fauvet et al. (2012): exponential radial profile, no halo.
- **Archimedean spiral, Stanev, TT**: zero field on the z-axis instead of NaN.
- **Derivatives**: `active_parameters` selects the columns, and columns follow `parameter_names` (columns were mislabelled or duplicated in some models); finite derivatives for SVT22 at z = 0 and JF12 on the z-axis.
- **Random fields**: samples realise exactly `rms(x)²` as local variance (amplitudes were arbitrary before); divergence cleaning and anisotropy work; `GaussianScalarField` uses `mu` and `sigma`; the ES profile follows hammurabiX. Random numbers are identical on all platforms, so samples for a given seed differ from earlier versions.

### Added

- Models: UF26 random field, XH24 halo, JF12 Planck variants (`Planck12b`, `Planck12c`, regular and random), Han XH24 disk variant.
- `PointCloud` grids for evaluation at arbitrary positions.
- `derivative` on grids and point clouds, shape `(3, ..., n_active)`.
- `interpolate`: linear or nearest-grid-point interpolation of data on a `RegularGrid` (e.g. a random sample) at arbitrary positions.
- `rms`, `variance` (and `mean` for scalar fields) for random fields.
- Eigen 5 support (Eigen 3.4 and 5.0 tested).
- Installable C++ library with a CMake package (`ImagineModels::ImagineModels`) and pkg-config file.

### API changes

- Grids are values passed to `evaluate(grid)` / `sample(grid, seed)`; `on_grid` and the grid constructors of the models are removed. Python results are NumPy arrays of shape `(3, nx, ny, nz)` or `(nx, ny, nz)`.
- Parameters: `parameter_names`, `parameters` (dict), `active_parameters`.
- Published variants are selected with `set_model` (UF24, TF17, Pshirkov, JF12, Han, UF26), which also loads their parameters.
- Classes renamed: `JF12RegularField` → `JF12MagneticField`, `SVT22` → `SVT22MagneticField` (Python, now equal to the C++ names), `UFMagneticField` → `UF24MagneticField`, `TFMagneticField` → `TF17MagneticField` (C++ and Python).
- Pshirkov: `set_model("ASS" | "BSS")` and `useDisk` replace `useASS` / `useBSS`.
- All angle inputs are in degrees. Changed from radians: UF24 `fDiskPhase1-3`, `fDiskPitch`, `fSpurCenter`, `fSpurLength`, `fSpurWidth`; Stanev `b_phi0`; YMW16 `t3_thmin`. Derivatives w.r.t. these parameters are per degree.
- YMW16: `t3_thmin`, `t3_tan_pitch`, `t3_cos_pitch` replace `t3_phimin` / `t3_tpitch`.
- Fauvet: parameters `b_b0, b_RB, b_Rsun, b_z0, b_r0, b_p, b_chi0`.
- ES random field: parameters `b0`, `observer`. `GaussianScalarField`: `mu`, `sigma`.
- Jaffe: `arm_num` outside 2–4 raises an error.
- C++: constants in `imagine::units` (`ImagineModels/units.h`; the `num`, `astro` and `cgs` namespaces are removed); all code in `namespace imagine`; headers included as `"ImagineModels/X.h"`; feature macros `IMAGINE_HAS_AUTODIFF`, `IMAGINE_HAS_FFTW`. Header files renamed to one name per model (e.g. `JF12.h`, `UF24.h`, `YMW16.h`).
- matplotlib is optional (extra `plot`) and no longer imported with the package.
