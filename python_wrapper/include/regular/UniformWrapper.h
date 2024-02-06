#ifndef UNIFORMWRAPPER_H
#define UNIFORMWRAPPER_H

#include <pybind11/pybind11.h>

#include "Uniform.h"

namespace py = pybind11;
using namespace pybind11::literals;

void Uniform(py::module_ &m)
{
    py::class_<UniformMagneticField, RegularVectorField>(m, "UniformMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("bx", &UniformMagneticField::bx)
        .def_readwrite("by", &UniformMagneticField::by)
        .def_readwrite("bz", &UniformMagneticField::bz)
#if autodiff_FOUND
        .def_readwrite("active_diff", &UniformMagneticField::active_diff)
        .def_readonly("all_diff", &UniformMagneticField::all_diff)

        .def("derivative", [](UniformMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif

        .def("at_position", [](UniformMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);

    py::class_<UniformDensityField, RegularScalarField>(m, "UniformDensityField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("n0", &UniformDensityField::n0)

#if autodiff_FOUND
        .def_readwrite("active_diff", &UniformDensityField::active_diff)
        .def_readonly("all_diff", &UniformDensityField::all_diff)

        .def("derivative", [](UniformDensityField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif

        .def("at_position", &UniformDensityField::at_position, "x"_a, "y"_a, "z"_a, py::return_value_policy::move);
}

#endif
