#ifndef REGULARFIELDBASES_H
#define REGULARFIELDBASES_H

#include <pybind11/pybind11.h>

#include "../regular_trampoline.h"
#include "../grid_bindings.h"

namespace py = pybind11;
using namespace pybind11::literals;

void RegularFieldBases(py::module_ &m) {
    py::class_<RegularVectorField, PyRegularVectorField>(m, "RegularVectorField")
        .def(py::init<>())
        .def("evaluate", [](const RegularVectorField &self, const RegularGrid &grid) { return to_numpy(self.evaluate(grid)); }, "grid"_a)
        .def("evaluate", [](const RegularVectorField &self, const IrregularGrid &grid) { return to_numpy(self.evaluate(grid)); }, "grid"_a);

    py::class_<RegularScalarField, PyRegularScalarField>(m, "RegularScalarField")
        .def(py::init<>())
        .def("evaluate", [](const RegularScalarField &self, const RegularGrid &grid) { return to_numpy(self.evaluate(grid)); }, "grid"_a)
        .def("evaluate", [](const RegularScalarField &self, const IrregularGrid &grid) { return to_numpy(self.evaluate(grid)); }, "grid"_a);
}

#endif
