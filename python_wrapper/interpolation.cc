#include <string>

#include "ImagineModels/Interpolation.h"
#include "bindings.h"

namespace {

using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

Interpolation to_method(const std::string &method) {
    if (method == "linear")
        return Interpolation::linear;
    if (method == "nearest")
        return Interpolation::nearest;
    throw GridException("interpolate: method must be 'linear' or 'nearest'.");
}

template <int N> GridData<N> to_grid_data(const DoubleArray &array) {
    std::array<int, 3> shape;
    for (int a = 0; a < 3; ++a)
        shape[a] = int(array.shape(array.ndim() - 3 + a));
    GridData<N> data(shape);
    std::copy(array.data(), array.data() + array.size(), data.data.begin());
    return data;
}

py::array_t<double> interpolate_points(const DoubleArray &data, const RegularGrid &grid, const PointCloud &points,
                                       const std::string &method, bool nan_outside) {
    if (data.ndim() == 4 && data.shape(0) == 3)
        return to_numpy_points(interpolate(to_grid_data<3>(data), grid, points, to_method(method), nan_outside));
    if (data.ndim() == 3)
        return to_numpy_points(interpolate(to_grid_data<1>(data), grid, points, to_method(method), nan_outside));
    throw GridException("interpolate: data must have shape (3, nx, ny, nz) or (nx, ny, nz).");
}

}

void bind_interpolation(py::module_ &m) {
    m.def("interpolate", &interpolate_points, "data"_a, "grid"_a, "points"_a, "method"_a = "linear",
          "nan_outside"_a = false, doc::interpolate);
    m.def(
        "interpolate",
        [](const DoubleArray &data, const RegularGrid &grid, const py::object &x, const py::object &y,
           const py::object &z, const std::string &method, bool nan_outside) {
            py::tuple broadcast = py::module_::import("numpy").attr("broadcast_arrays")(x, y, z);
            DoubleArray xs = broadcast[0].cast<DoubleArray>();
            DoubleArray ys = broadcast[1].cast<DoubleArray>();
            DoubleArray zs = broadcast[2].cast<DoubleArray>();
            const PointCloud points(std::vector<double>(xs.data(), xs.data() + xs.size()),
                                    std::vector<double>(ys.data(), ys.data() + ys.size()),
                                    std::vector<double>(zs.data(), zs.data() + zs.size()));
            py::array_t<double> flat = interpolate_points(data, grid, points, method, nan_outside);
            std::vector<py::ssize_t> shape(xs.shape(), xs.shape() + xs.ndim());
            if (flat.ndim() == 2)
                shape.insert(shape.begin(), 3);
            return flat.reshape(shape);
        },
        "data"_a, "grid"_a, "x"_a, "y"_a, "z"_a, "method"_a = "linear", "nan_outside"_a = false);
}
