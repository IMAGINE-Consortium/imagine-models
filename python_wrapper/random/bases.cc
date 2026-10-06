#include "../bindings.h"
#include "../random_trampoline.h"

template <typename Field, typename PyClass>
void bind_statistics(PyClass &cls)
{
    cls.def("rms", [](const Field &self, const RegularGrid &grid) { return to_numpy(self.evaluate_rms(grid)); }, "grid"_a)
        .def("rms", [](const Field &self, const IrregularGrid &grid) { return to_numpy(self.evaluate_rms(grid)); }, "grid"_a)
        .def("rms", [](const Field &self, const py::object &x, const py::object &y, const py::object &z) {
            return map_positions([&](double a, double b, double c) { return self.rms(a, b, c); }, x, y, z); }, "x"_a, "y"_a, "z"_a)
        .def("variance", [](const Field &self, const py::object &x, const py::object &y, const py::object &z) {
            return map_positions([&](double a, double b, double c) { return self.variance(a, b, c); }, x, y, z); }, "x"_a, "y"_a, "z"_a)
        .def("spectrum", &Field::spectrum, "abs_k"_a)
        .def_readwrite("apply_spectrum", &Field::apply_spectrum);
}

void bind_random_bases(py::module_ &m)
{
    py::class_<RandomVectorField, PyRandomVectorField> vector(m, "RandomVectorField");
    vector.def(py::init<>())
        .def("sample", [](const RandomVectorField &self, const RegularGrid &grid, int seed) { return to_numpy(self.sample(grid, seed)); }, "grid"_a, "seed"_a)
        .def("random_numbers", [](const RandomVectorField &self, const RegularGrid &grid, int seed) { return to_numpy(self.random_numbers(grid, seed)); }, "grid"_a, "seed"_a)
        .def("anisotropy_direction", &RandomVectorField::anisotropy_direction, "x"_a, "y"_a, "z"_a)
        .def_readwrite("clean_divergence", &RandomVectorField::clean_divergence)
        .def_readwrite("apply_anisotropy", &RandomVectorField::apply_anisotropy)
        .def_readwrite("anisotropy_rho", &RandomVectorField::anisotropy_rho);
    bind_statistics<RandomVectorField>(vector);

    py::class_<RandomScalarField, PyRandomScalarField> scalar(m, "RandomScalarField");
    scalar.def(py::init<>())
        .def("sample", [](const RandomScalarField &self, const RegularGrid &grid, int seed) { return to_numpy(self.sample(grid, seed)); }, "grid"_a, "seed"_a)
        .def("random_numbers", [](const RandomScalarField &self, const RegularGrid &grid, int seed) { return to_numpy(self.random_numbers(grid, seed)); }, "grid"_a, "seed"_a)
        .def("mean", [](const RandomScalarField &self, const py::object &x, const py::object &y, const py::object &z) {
            return map_positions([&](double a, double b, double c) { return self.mean(a, b, c); }, x, y, z); }, "x"_a, "y"_a, "z"_a);
    bind_statistics<RandomScalarField>(scalar);
}
