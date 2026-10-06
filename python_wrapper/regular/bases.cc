#include "../bindings.h"
#include "../regular_trampoline.h"

void bind_regular_bases(py::module_ &m)
{
    py::class_<RegularVectorField, PyRegularVectorField>(m, "RegularVectorField")
        .def(py::init<>())
        .def("evaluate", [](const RegularVectorField &self, const RegularGrid &grid) { return to_numpy(self.evaluate(grid)); }, "grid"_a)
        .def("evaluate", [](const RegularVectorField &self, const IrregularGrid &grid) { return to_numpy(self.evaluate(grid)); }, "grid"_a);

    py::class_<RegularScalarField, PyRegularScalarField>(m, "RegularScalarField")
        .def(py::init<>())
        .def("evaluate", [](const RegularScalarField &self, const RegularGrid &grid) { return to_numpy(self.evaluate(grid)); }, "grid"_a)
        .def("evaluate", [](const RegularScalarField &self, const IrregularGrid &grid) { return to_numpy(self.evaluate(grid)); }, "grid"_a);
}
