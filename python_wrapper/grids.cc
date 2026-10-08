#include "bindings.h"

void bind_grids(py::module_ &m) {
    py::class_<RegularGrid>(m, "RegularGrid", "Regular grid: shape, reference_point (first point), increment.")
        .def(py::init<const std::array<int, 3> &, const std::array<double, 3> &, const std::array<double, 3> &>(),
             py::arg("shape").noconvert(), "reference_point"_a, "increment"_a)
        .def_readonly("shape", &RegularGrid::shape)
        .def_readonly("reference_point", &RegularGrid::reference_point)
        .def_readonly("increment", &RegularGrid::increment)
        .def_property_readonly("size", &RegularGrid::size)
        .def("__repr__", [](const RegularGrid &g) {
            return py::str("RegularGrid(shape={}, reference_point={}, increment={})")
                .format(g.shape, g.reference_point, g.increment);
        });

    py::class_<IrregularGrid>(m, "IrregularGrid", "Grid spanned by arbitrary x, y and z axes.")
        .def(py::init<const std::vector<double> &, const std::vector<double> &, const std::vector<double> &>(), "x"_a,
             "y"_a, "z"_a)
        .def_property_readonly("x", [](const IrregularGrid &g) { return py::array_t<double>(g.x.size(), g.x.data()); })
        .def_property_readonly("y", [](const IrregularGrid &g) { return py::array_t<double>(g.y.size(), g.y.data()); })
        .def_property_readonly("z", [](const IrregularGrid &g) { return py::array_t<double>(g.z.size(), g.z.data()); })
        .def_property_readonly("shape", &IrregularGrid::shape)
        .def_property_readonly("size", &IrregularGrid::size)
        .def("__repr__", [](const IrregularGrid &g) { return py::str("IrregularGrid(shape={})").format(g.shape()); });

    py::class_<PointCloud>(m, "PointCloud", "Unstructured set of N positions.")
        .def(py::init<const std::vector<double> &, const std::vector<double> &, const std::vector<double> &>(), "x"_a,
             "y"_a, "z"_a)
        .def_static(
            "from_positions",
            [](const py::array_t<double, py::array::c_style | py::array::forcecast> &positions) {
                if (positions.ndim() != 2 || positions.shape(1) != 3)
                    throw GridException("PointCloud.from_positions: positions must have shape (N, 3).");
                const py::ssize_t n = positions.shape(0);
                std::vector<double> x(n), y(n), z(n);
                auto p = positions.unchecked<2>();
                for (py::ssize_t i = 0; i < n; ++i) {
                    x[i] = p(i, 0);
                    y[i] = p(i, 1);
                    z[i] = p(i, 2);
                }
                return PointCloud(x, y, z);
            },
            "positions"_a)
        .def_property_readonly("x", [](const PointCloud &g) { return py::array_t<double>(g.x.size(), g.x.data()); })
        .def_property_readonly("y", [](const PointCloud &g) { return py::array_t<double>(g.y.size(), g.y.data()); })
        .def_property_readonly("z", [](const PointCloud &g) { return py::array_t<double>(g.z.size(), g.z.data()); })
        .def_property_readonly("size", &PointCloud::size)
        .def("__len__", &PointCloud::size)
        .def("__repr__", [](const PointCloud &g) { return py::str("PointCloud(size={})").format(g.size()); });

    py::register_exception<GridException>(m, "GridError", PyExc_ValueError);
}
