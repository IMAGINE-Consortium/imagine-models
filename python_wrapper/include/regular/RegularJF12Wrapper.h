#ifndef REGULARJF12WRAPPER_H
#define REGULARJF12WRAPPER_H

#include <pybind11/pybind11.h>

#include "RegularJF12.h"

void RegularJF12(py::module_ &m)
{
    py::class_<JF12MagneticField, RegularVectorField>(m, "JF12RegularField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("b_arm_1", &JF12MagneticField::b_arm_1)
        .def_readwrite("b_arm_2", &JF12MagneticField::b_arm_2)
        .def_readwrite("b_arm_3", &JF12MagneticField::b_arm_3)
        .def_readwrite("b_arm_4", &JF12MagneticField::b_arm_4)
        .def_readwrite("b_arm_5", &JF12MagneticField::b_arm_5)
        .def_readwrite("b_arm_6", &JF12MagneticField::b_arm_6)
        .def_readwrite("b_arm_7", &JF12MagneticField::b_arm_7)
        .def_readwrite("b_ring", &JF12MagneticField::b_ring)
        .def_readwrite("h_disk", &JF12MagneticField::h_disk)
        .def_readwrite("w_disk", &JF12MagneticField::w_disk)
        .def_readwrite("do_halo", &JF12MagneticField::do_halo)
        .def_readwrite("Bn", &JF12MagneticField::Bn)
        .def_readwrite("Bs", &JF12MagneticField::Bs)
        .def_readwrite("rn", &JF12MagneticField::rn)
        .def_readwrite("rs", &JF12MagneticField::rs)
        .def_readwrite("wh", &JF12MagneticField::wh)
        .def_readwrite("z0", &JF12MagneticField::z0)
        .def_readwrite("do_X", &JF12MagneticField::do_X)
        .def_readwrite("B0_X", &JF12MagneticField::B0_X)
        .def_readwrite("Xtheta_const", &JF12MagneticField::Xtheta_const)
        .def_readwrite("rpc_X", &JF12MagneticField::rpc_X)
        .def_readwrite("r0_X", &JF12MagneticField::r0_X)
#if autodiff_FOUND
        .def_readwrite("active_diff", &JF12MagneticField::active_diff)
        .def_readonly("all_diff", &JF12MagneticField::all_diff)

        .def("derivative", [](JF12MagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif

        .def("at_position", [](JF12MagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif