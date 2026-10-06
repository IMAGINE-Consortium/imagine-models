#ifndef HELIXWRAPPER_H
#define HELIXWRAPPER_H

#include <pybind11/pybind11.h>

#include "ImagineModels/Helix.h"

namespace py = pybind11;
using namespace pybind11::literals;

void Helix(py::module_ &m)
{

    py::class_<HelixMagneticField, RegularVectorField>(m, "HelixMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("ampx", &HelixMagneticField::ampx)
        .def_readwrite("ampy", &HelixMagneticField::ampy)
        .def_readwrite("ampz", &HelixMagneticField::ampz)
        .def_readwrite("rmin", &HelixMagneticField::rmin)
        .def_readwrite("rmax", &HelixMagneticField::rmax)
#if IMAGINE_HAS_AUTODIFF
        .def_readwrite("active_diff", &HelixMagneticField::active_diff)
        .def_readonly("all_diff", &HelixMagneticField::all_diff)

        .def("derivative", [](HelixMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif

    .def("at_position", [](HelixMagneticField &self, double x, double y, double z)
        {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
        },
        "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif