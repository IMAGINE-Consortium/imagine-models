#ifndef GRID_BINDINGS_H
#define GRID_BINDINGS_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "ImagineModels/Grid.h"

namespace py = pybind11;
using namespace pybind11::literals;

template <int N>
py::array_t<double> to_numpy(GridData<N> &&grid_data) {
  auto owned = new GridData<N>(std::move(grid_data));
  py::capsule owner(owned, [](void *p) { delete static_cast<GridData<N> *>(p); });
  std::vector<py::ssize_t> shape;
  if (N > 1)
    shape.push_back(N);
  for (int s : owned->shape)
    shape.push_back(s);
  return py::array_t<double>(shape, owned->data.data(), owner);
}

void Grids(py::module_ &m) {
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

#endif
