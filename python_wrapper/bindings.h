#pragma once

#include <utility>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ImagineModels/Grid.h"
#include "ImagineModels/config.h"

#if IMAGINE_HAS_AUTODIFF
#include <pybind11/eigen.h>
#endif

namespace py = pybind11;
using namespace pybind11::literals;
using namespace imagine;

namespace doc {
inline constexpr const char *evaluate = "Field on a RegularGrid, IrregularGrid or PointCloud.";
inline constexpr const char *at_position = "Field at one position (kpc).";
inline constexpr const char *at_positions = "Field at many positions, with NumPy broadcasting.";
inline constexpr const char *derivative =
    "Jacobian w.r.t. the active parameters at one position or on a grid (parameters on the last axis).";
inline constexpr const char *active_parameters = "Parameters included in derivative, in this order.";
inline constexpr const char *parameters = "All parameters as a dict; assigning updates the given ones.";
inline constexpr const char *parameter_names = "Names of the model parameters, in order.";
inline constexpr const char *set_model = "Select a published variant and load its parameters.";
inline constexpr const char *sample = "Random realisation on a RegularGrid; use interpolate for other positions.";
inline constexpr const char *interpolate =
    "Linear or nearest-grid-point interpolation of grid data at positions; raises outside the grid unless nan_outside.";
inline constexpr const char *random_numbers = "Unit-variance Gaussian random field before scaling.";
inline constexpr const char *rms = "Expected rms amplitude at positions or on a grid.";
inline constexpr const char *variance = "Expected variance at positions.";
inline constexpr const char *mean = "Expected mean at positions.";
inline constexpr const char *spectrum = "Power spectrum at wave number |k|.";
}

template <int N> py::array_t<double> to_numpy(GridData<N> &&grid_data) {
    auto owned = new GridData<N>(std::move(grid_data));
    py::capsule owner(owned, [](void *p) { delete static_cast<GridData<N> *>(p); });
    std::vector<py::ssize_t> shape;
    if (N > 1)
        shape.push_back(N);
    for (int s : owned->shape)
        shape.push_back(s);
    return py::array_t<double>(shape, owned->data.data(), owner);
}

template <int N> py::array_t<double> to_numpy_points(GridData<N> &&grid_data) {
    auto owned = new GridData<N>(std::move(grid_data));
    py::capsule owner(owned, [](void *p) { delete static_cast<GridData<N> *>(p); });
    std::vector<py::ssize_t> shape;
    if (N > 1)
        shape.push_back(N);
    shape.push_back(owned->shape[0]);
    return py::array_t<double>(shape, owned->data.data(), owner);
}

template <int N> py::array_t<double> to_numpy_jacobian(JacobianData<N> &&jacobian, bool points) {
    auto owned = new JacobianData<N>(std::move(jacobian));
    py::capsule owner(owned, [](void *p) { delete static_cast<JacobianData<N> *>(p); });
    std::vector<py::ssize_t> shape;
    if (N > 1)
        shape.push_back(N);
    for (int a = 0; a < (points ? 1 : 3); ++a)
        shape.push_back(owned->shape[a]);
    shape.push_back(py::ssize_t(owned->columns));
    return py::array_t<double>(shape, owned->data.data(), owner);
}

template <typename F>
py::array_t<double> map_positions(F &&f, const py::object &x, const py::object &y, const py::object &z) {
    using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;
    py::tuple broadcast = py::module_::import("numpy").attr("broadcast_arrays")(x, y, z);
    DoubleArray xs = broadcast[0].cast<DoubleArray>();
    DoubleArray ys = broadcast[1].cast<DoubleArray>();
    DoubleArray zs = broadcast[2].cast<DoubleArray>();
    py::array_t<double> out(std::vector<py::ssize_t>(xs.shape(), xs.shape() + xs.ndim()));
    double *o = out.mutable_data();
    for (py::ssize_t i = 0; i < xs.size(); ++i)
        o[i] = f(xs.data()[i], ys.data()[i], zs.data()[i]);
    return out;
}

void bind_grids(py::module_ &m);
void bind_interpolation(py::module_ &m);
void bind_regular_bases(py::module_ &m);
void bind_archimedes(py::module_ &m);
void bind_fauvet(py::module_ &m);
void bind_han(py::module_ &m);
void bind_hmr(py::module_ &m);
void bind_helix(py::module_ &m);
void bind_jaffe(py::module_ &m);
void bind_kst24(py::module_ &m);
void bind_ne2025(py::module_ &m);
void bind_plane_parallel(py::module_ &m);
void bind_yt20(py::module_ &m);
void bind_pshirkov(py::module_ &m);
void bind_jf12(py::module_ &m);
void bind_stanev(py::module_ &m);
void bind_sun(py::module_ &m);
void bind_svt22(py::module_ &m);
void bind_tf17(py::module_ &m);
void bind_tt(py::module_ &m);
void bind_uf24(py::module_ &m);
void bind_uniform(py::module_ &m);
void bind_wmap(py::module_ &m);
void bind_ymw16(py::module_ &m);
void bind_xh24(py::module_ &m);

#if IMAGINE_HAS_FFTW
void bind_random_bases(py::module_ &m);
void bind_es_random(py::module_ &m);
void bind_gaussian_scalar(py::module_ &m);
void bind_lognormal(py::module_ &m);
void bind_jf12_random(py::module_ &m);
void bind_uf26_random(py::module_ &m);
#endif
