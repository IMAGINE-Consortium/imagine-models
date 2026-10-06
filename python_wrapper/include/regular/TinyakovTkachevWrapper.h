#ifndef TTWRAPPER_H
#define TTWRAPPER_H

#include <pybind11/pybind11.h>

#include "ImagineModels/TinyakovTkachev.h"

namespace py = pybind11;
using namespace pybind11::literals;

void TinyakovTkachev(py::module_ &m)
{
    py::class_<TTMagneticField, RegularVectorField>(m, "TTMagneticField")
        .def(py::init<>())
        .def_readwrite("b_b0", &TTMagneticField::b_b0)
        .def_readwrite("b_Rsun", &TTMagneticField::b_Rsun)
        .def_readwrite("b_r_max", &TTMagneticField::b_r_max)
        .def_readwrite("b_r_min", &TTMagneticField::b_r_min)
        .def_readwrite("b_d", &TTMagneticField::b_d)
        .def_readwrite("b_z0", &TTMagneticField::b_z0)
        .def_readwrite("b_p", &TTMagneticField::b_p)
#if IMAGINE_HAS_AUTODIFF
        .def_readwrite("active_diff", &TTMagneticField::active_diff)
        .def_readonly("all_diff", &TTMagneticField::all_diff)

        .def("derivative", [](TTMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif
        .def(
            "at_position", [](TTMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            }, 
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif