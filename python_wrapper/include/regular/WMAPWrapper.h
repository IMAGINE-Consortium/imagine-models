#ifndef WMAPWRAPPER_H
#define WMAPWRAPPER_H

#include <pybind11/pybind11.h>

#include "ImagineModels/WMAP.h"

namespace py = pybind11;
using namespace pybind11::literals;

void WMAP(py::module_ &m) {
    py::class_<WMAPMagneticField, RegularVectorField>(m, "WMAPMagneticField")
        .def(py::init<>())

        .def_readwrite("b_b0", &WMAPMagneticField::b_b0)
        .def_readwrite("b_r0", &WMAPMagneticField::b_r0)
        .def_readwrite("b_r_max", &WMAPMagneticField::b_r_max)
        .def_readwrite("b_r_min", &WMAPMagneticField::b_r_min)
        .def_readwrite("b_Rsun", &WMAPMagneticField::b_Rsun)
        .def_readwrite("b_z0", &WMAPMagneticField::b_z0)
        
        .def_readwrite("b_psi0", &WMAPMagneticField::b_psi0)
        .def_readwrite("b_psi1", &WMAPMagneticField::b_psi1)
        .def_readwrite("b_xsi0", &WMAPMagneticField::b_xsi0)
        .def_readwrite("b_anti", &WMAPMagneticField::anti)
#if IMAGINE_HAS_AUTODIFF
        .def_readwrite("active_diff", &WMAPMagneticField::active_diff)
        .def_readonly("all_diff", &WMAPMagneticField::all_diff)

        .def("derivative", [](WMAPMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif
        .def("at_position", [](WMAPMagneticField &self, double x, double y, double z)  
            {
                vector f = self.at_position(x, y, z);
                return std::make_tuple(f[0], f[1], f[2]);
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif