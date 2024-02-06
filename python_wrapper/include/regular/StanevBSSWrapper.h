#ifndef STANEVBSSWRAPPER_H
#define STANEVBSSWRAPPER_H

#include <pybind11/pybind11.h>

#include "StanevBSS.h"

namespace py = pybind11;
using namespace pybind11::literals;

void StanevBSS(py::module_ &m)
{
    py::class_<StanevBSSMagneticField, RegularVectorField>(m, "StanevBSSMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("b_r0", &StanevBSSMagneticField::b_r0)
        .def_readwrite("b_phi0", &StanevBSSMagneticField::b_phi0)
        .def_readwrite("b_r_max", &StanevBSSMagneticField::b_r_max)
        .def_readwrite("b_r_min", &StanevBSSMagneticField::b_r_min)
        .def_readwrite("b_Rsun", &StanevBSSMagneticField::b_Rsun)
        .def_readwrite("b_z01", &StanevBSSMagneticField::b_z01)
        .def_readwrite("b_z02", &StanevBSSMagneticField::b_z02)
        .def_readwrite("b_z0_border", &StanevBSSMagneticField::b_z0_border)
        .def_readwrite("b_p", &StanevBSSMagneticField::b_p)
#if autodiff_FOUND
        .def_readwrite("active_diff", &StanevBSSMagneticField::active_diff)
        .def_readonly("all_diff", &StanevBSSMagneticField::all_diff)

        .def("derivative", [](StanevBSSMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif

        .def("at_position", [](StanevBSSMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif