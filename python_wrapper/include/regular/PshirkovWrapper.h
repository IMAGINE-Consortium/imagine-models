#ifndef PSHIRKOVWRAPPER_H
#define PSHIRKOVWRAPPER_H

#include <pybind11/pybind11.h>

#include "Pshirkov.h"

void Pshirkov(py::module_ &m)
{
    py::class_<PshirkovMagneticField, RegularVectorField>(m, "PshirkovMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("useASS", &PshirkovMagneticField::useASS)
        .def_readwrite("useBSS", &PshirkovMagneticField::useBSS)
        .def_readwrite("useHalo", &PshirkovMagneticField::useHalo)

        .def_readwrite("pitch", &PshirkovMagneticField::pitch)
        .def_readwrite("d", &PshirkovMagneticField::d)
        .def_readwrite("R_sun", &PshirkovMagneticField::R_sun)
        .def_readwrite("R_c", &PshirkovMagneticField::R_c)
        .def_readwrite("z0_D", &PshirkovMagneticField::z0_D)
        .def_readwrite("B0_D", &PshirkovMagneticField::B0_D)
        .def_readwrite("z0_H", &PshirkovMagneticField::z0_H)
        .def_readwrite("R0_H", &PshirkovMagneticField::R0_H)
        .def_readwrite("B0_Hn", &PshirkovMagneticField::B0_Hn)
        .def_readwrite("B0_Hs", &PshirkovMagneticField::B0_Hs)
        .def_readwrite("z11_H", &PshirkovMagneticField::z11_H)
        .def_readwrite("z12_H", &PshirkovMagneticField::z12_H)

#if autodiff_FOUND
        .def_readwrite("active_diff", &PshirkovMagneticField::active_diff)
        .def_readonly("all_diff", &PshirkovMagneticField::all_diff)

        .def("derivative", [](PshirkovMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif
        .def(
            "at_position", [](PshirkovMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif