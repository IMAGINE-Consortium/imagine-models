#include "bindings.h"

void bind_grids(py::module_ &m)
{
  py::class_<RegularGrid>(m, "RegularGrid")
      .def(py::init<const std::array<int, 3> &, const std::array<double, 3> &, const std::array<double, 3> &>(),
           py::arg("shape").noconvert(), "reference_point"_a, "increment"_a)
      .def_readonly("shape", &RegularGrid::shape)
      .def_readonly("reference_point", &RegularGrid::reference_point)
      .def_readonly("increment", &RegularGrid::increment)
      .def_property_readonly("size", &RegularGrid::size)
      .def("__repr__", [](const RegularGrid &g) {
        return py::str("RegularGrid(shape={}, reference_point={}, increment={})").format(g.shape, g.reference_point, g.increment);
      });

  py::class_<IrregularGrid>(m, "IrregularGrid")
      .def(py::init<const std::vector<double> &, const std::vector<double> &, const std::vector<double> &>(), "x"_a, "y"_a, "z"_a)
      .def_property_readonly("x", [](const IrregularGrid &g) { return py::array_t<double>(g.x.size(), g.x.data()); })
      .def_property_readonly("y", [](const IrregularGrid &g) { return py::array_t<double>(g.y.size(), g.y.data()); })
      .def_property_readonly("z", [](const IrregularGrid &g) { return py::array_t<double>(g.z.size(), g.z.data()); })
      .def_property_readonly("shape", &IrregularGrid::shape)
      .def_property_readonly("size", &IrregularGrid::size)
      .def("__repr__", [](const IrregularGrid &g) {
        return py::str("IrregularGrid(shape={})").format(g.shape());
      });

  py::register_exception<GridException>(m, "GridError", PyExc_ValueError);
}
