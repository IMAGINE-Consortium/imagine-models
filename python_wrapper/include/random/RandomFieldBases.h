#ifndef RANDOMFIELDBASES_H
#define RANDOMFIELDBASES_H

#include <pybind11/pybind11.h>

#include "../random_trampoline.h"
#include "../grid_bindings.h"

namespace py = pybind11;
using namespace pybind11::literals;

void RandomFieldBases(py::module_ &m) {
    py::class_<RandomVectorField, PyRandomVectorField>(m, "RandomVectorField")
        .def(py::init<>())
        .def("sample", [](const RandomVectorField &self, const RegularGrid &grid, int seed) { return to_numpy(self.sample(grid, seed)); }, "grid"_a, "seed"_a)
        .def("random_numbers", [](const RandomVectorField &self, const RegularGrid &grid, int seed) { return to_numpy(self.random_numbers(grid, seed)); }, "grid"_a, "seed"_a)
        .def("profile", [](const RandomVectorField &self, const RegularGrid &grid) { return to_numpy(self.profile(grid)); }, "grid"_a)
        .def("profile", [](const RandomVectorField &self, const IrregularGrid &grid) { return to_numpy(self.profile(grid)); }, "grid"_a);

    py::class_<RandomScalarField, PyRandomScalarField>(m, "RandomScalarField")
        .def(py::init<>())
        .def("sample", [](const RandomScalarField &self, const RegularGrid &grid, int seed) { return to_numpy(self.sample(grid, seed)); }, "grid"_a, "seed"_a)
        .def("random_numbers", [](const RandomScalarField &self, const RegularGrid &grid, int seed) { return to_numpy(self.random_numbers(grid, seed)); }, "grid"_a, "seed"_a)
        .def("profile", [](const RandomScalarField &self, const RegularGrid &grid) { return to_numpy(self.profile(grid)); }, "grid"_a)
        .def("profile", [](const RandomScalarField &self, const IrregularGrid &grid) { return to_numpy(self.profile(grid)); }, "grid"_a);
}

#endif
