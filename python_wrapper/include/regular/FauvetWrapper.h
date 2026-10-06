#ifndef FAUVETWRAPPER_H
#define FAUVETWRAPPER_H

#include <pybind11/pybind11.h>

#include "ImagineModels/Fauvet.h"

namespace py = pybind11;
using namespace pybind11::literals;

void Fauvet(py::module_ &m)
{
    py::class_<FauvetMagneticField, RegularVectorField>(m, "FauvetMagneticField")
        .def(py::init<>())

        .def_readwrite("b_b0", &FauvetMagneticField::b_b0)
        .def_readwrite("b_r0", &FauvetMagneticField::b_r0)
        .def_readwrite("b_r_max", &FauvetMagneticField::b_r_max)
        .def_readwrite("b_r_min", &FauvetMagneticField::b_r_min)
        .def_readwrite("b_z0", &FauvetMagneticField::b_z0)
        .def_readwrite("b_p", &FauvetMagneticField::b_p)
        .def_readwrite("b_chi0", &FauvetMagneticField::b_chi0)

        .def_readwrite("h_b0", &FauvetMagneticField::h_b0)
        .def_readwrite("h_r0", &FauvetMagneticField::h_r0)
        .def_readwrite("h_z0", &FauvetMagneticField::h_z0)
        .def_readwrite("h_z1a", &FauvetMagneticField::h_z1a)
        .def_readwrite("h_z1b", &FauvetMagneticField::h_z1b)
#if IMAGINE_HAS_AUTODIFF
        .def_readwrite("active_diff", &FauvetMagneticField::active_diff)
        .def_readonly("all_diff", &FauvetMagneticField::all_diff)

        .def("derivative", [](FauvetMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif

        .def("at_position", [](FauvetMagneticField &self, double x, double y, double z)
            {
                vector f = self.at_position(x, y, z);
                return std::make_tuple(f[0], f[1], f[2]);
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif
