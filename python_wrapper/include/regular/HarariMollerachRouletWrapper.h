#ifndef HMRWRAPPER_H
#define HMRWRAPPER_H

#include <pybind11/pybind11.h>

#include "ImagineModels/HarariMollerachRoulet.h"

namespace py = pybind11;
using namespace pybind11::literals;

void HarariMollerachRoulet(py::module_ &m)
{
    py::class_<HMRMagneticField, RegularVectorField>(m, "HMRMagneticField")
        .def(py::init<>())

        .def_readwrite("b_Rsun", &HMRMagneticField::b_Rsun)
        .def_readwrite("b_r_max", &HMRMagneticField::b_r_max)
        .def_readwrite("b_z1", &HMRMagneticField::b_z1)
        .def_readwrite("b_z2", &HMRMagneticField::b_z2)
        .def_readwrite("b_r1", &HMRMagneticField::b_r1)
        .def_readwrite("b_p", &HMRMagneticField::b_p)
        .def_readwrite("b_epsilon0", &HMRMagneticField::b_epsilon0)
#if IMAGINE_HAS_AUTODIFF
        .def_readwrite("active_diff", &HMRMagneticField::active_diff)
        .def_readonly("all_diff", &HMRMagneticField::all_diff)

        .def("derivative", [](HMRMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif

        .def("at_position", [](HMRMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif