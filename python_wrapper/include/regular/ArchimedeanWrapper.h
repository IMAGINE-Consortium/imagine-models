#ifndef ARCHIWRAPPER_H
#define ARCHIWRAPPER_H

#include <pybind11/pybind11.h>

#include "Archimedes.h"

void Archimedes(py::module_ &m)
{
    py::class_<ArchimedeanMagneticField, RegularVectorField>(m, "ArchimedeanMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("R_0", &ArchimedeanMagneticField::R_0)
        .def_readwrite("Omega", &ArchimedeanMagneticField::Omega)
        .def_readwrite("v_w", &ArchimedeanMagneticField::v_w)
        .def_readwrite("B_0", &ArchimedeanMagneticField::B_0)

#if autodiff_FOUND
        .def_readwrite("active_diff", &ArchimedeanMagneticField::active_diff)
        .def_readonly("all_diff", &ArchimedeanMagneticField::all_diff)

        .def("derivative", [](ArchimedeanMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif
        .def(
            "at_position", [](ArchimedeanMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif