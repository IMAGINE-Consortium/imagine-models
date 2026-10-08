#include <type_traits>

#include "../bindings.h"
#include "../regular_trampoline.h"

using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

template <typename Field>
py::array_t<double> at_positions(const Field &self, const py::object &x, const py::object &y, const py::object &z) {
    constexpr bool is_vector = std::is_base_of_v<RegularVectorField, Field>;
    py::tuple broadcast = py::module_::import("numpy").attr("broadcast_arrays")(x, y, z);
    DoubleArray xs = broadcast[0].cast<DoubleArray>();
    DoubleArray ys = broadcast[1].cast<DoubleArray>();
    DoubleArray zs = broadcast[2].cast<DoubleArray>();

    std::vector<py::ssize_t> shape(xs.shape(), xs.shape() + xs.ndim());
    const py::ssize_t n = xs.size();
    if (is_vector)
        shape.insert(shape.begin(), 3);
    py::array_t<double> out(shape);
    double *o = out.mutable_data();
    const double *px = xs.data(), *py_ = ys.data(), *pz = zs.data();
    for (py::ssize_t i = 0; i < n; ++i) {
        if constexpr (is_vector) {
            Vec3<double> b = self.at_position(px[i], py_[i], pz[i]);
            o[i] = b[0];
            o[n + i] = b[1];
            o[2 * n + i] = b[2];
        } else {
            o[i] = self.at_position(px[i], py_[i], pz[i]);
        }
    }
    return out;
}

void bind_regular_bases(py::module_ &m) {
    py::class_<RegularVectorField, PyRegularVectorField>(m, "RegularVectorField")
        .def(py::init<>())
        .def(
            "evaluate",
            [](const RegularVectorField &self, const RegularGrid &grid) { return to_numpy(self.evaluate(grid)); },
            "grid"_a, doc::evaluate)
        .def(
            "evaluate",
            [](const RegularVectorField &self, const IrregularGrid &grid) { return to_numpy(self.evaluate(grid)); },
            "grid"_a)
        .def(
            "evaluate",
            [](const RegularVectorField &self, const PointCloud &grid) { return to_numpy_points(self.evaluate(grid)); },
            "grid"_a)
        .def("at_positions", &at_positions<RegularVectorField>, "x"_a, "y"_a, "z"_a, doc::at_positions);

    py::class_<RegularScalarField, PyRegularScalarField>(m, "RegularScalarField")
        .def(py::init<>())
        .def(
            "evaluate",
            [](const RegularScalarField &self, const RegularGrid &grid) { return to_numpy(self.evaluate(grid)); },
            "grid"_a, doc::evaluate)
        .def(
            "evaluate",
            [](const RegularScalarField &self, const IrregularGrid &grid) { return to_numpy(self.evaluate(grid)); },
            "grid"_a)
        .def(
            "evaluate",
            [](const RegularScalarField &self, const PointCloud &grid) { return to_numpy_points(self.evaluate(grid)); },
            "grid"_a)
        .def("at_positions", &at_positions<RegularScalarField>, "x"_a, "y"_a, "z"_a, doc::at_positions);
}
