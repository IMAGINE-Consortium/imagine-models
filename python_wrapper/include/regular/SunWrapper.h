#ifndef SUNWRAPPER_H
#define SUNWRAPPER_H

#include <pybind11/pybind11.h>

#include "ImagineModels/Sun.h"

void Sun2008(py::module_ &m)
{
    py::class_<SunMagneticField, RegularVectorField>(m, "SunMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("b_B0", &SunMagneticField::b_B0)
        .def_readwrite("b_Rsun", &SunMagneticField::b_Rsun)
        .def_readwrite("b_R0", &SunMagneticField::b_R0)
        .def_readwrite("b_z0", &SunMagneticField::b_z0)
        .def_readwrite("b_Rc", &SunMagneticField::b_Rc)
        .def_readwrite("b_Bc", &SunMagneticField::b_Bc)
        .def_readwrite("b_p", &SunMagneticField::b_p)

        .def_readwrite("bH_B0", &SunMagneticField::bH_B0)
        .def_readwrite("bH_z0", &SunMagneticField::bH_z0)
        .def_readwrite("bH_z1a", &SunMagneticField::bH_z1a)
        .def_readwrite("bH_z1b", &SunMagneticField::bH_z1b)
        .def_readwrite("bH_R0", &SunMagneticField::bH_R0)
#if IMAGINE_HAS_AUTODIFF
        .def_readwrite("active_diff", &SunMagneticField::active_diff)
        .def_readonly("all_diff", &SunMagneticField::all_diff)

        .def("derivative", [](SunMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif
        .def(
            "at_position", [](SunMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif