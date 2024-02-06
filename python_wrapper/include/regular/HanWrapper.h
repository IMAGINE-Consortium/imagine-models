#ifndef HANWRAPPER_H
#define HANWRAPPER_H

#include <pybind11/pybind11.h>

#include "Han.h"

void Han2018(py::module_ &m)
{
    py::class_<HanMagneticField, RegularVectorField>(m, "HanMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("B_p", &HanMagneticField::B_p)
        .def_readwrite("A", &HanMagneticField::A)
        .def_readwrite("H", &HanMagneticField::H)
        .def_readwrite("B_s1", &HanMagneticField::B_s1)
        .def_readwrite("B_s2", &HanMagneticField::B_s2)
        .def_readwrite("B_s3", &HanMagneticField::B_s3)
        .def_readwrite("B_s4", &HanMagneticField::B_s4)
        .def_readwrite("B_s5", &HanMagneticField::B_s5)
        .def_readwrite("B_s6", &HanMagneticField::B_s6)

#if autodiff_FOUND
        .def_readwrite("active_diff", &HanMagneticField::active_diff)
        .def_readonly("all_diff", &HanMagneticField::all_diff)

        .def("derivative", [](HanMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif
        .def(
            "at_position", [](HanMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif