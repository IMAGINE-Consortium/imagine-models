#ifndef SVT22WRAPPER_H
#define SVT22WRAPPER_H

#include <pybind11/pybind11.h>

#include "SVT22.h"

void SVT22(py::module_ &m)
{
    py::class_<SVT22MagneticField, RegularVectorField>(m, "SVT22")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("B_val", &SVT22MagneticField::B_val)
        .def_readwrite("r_cut", &SVT22MagneticField::r_cut)
        .def_readwrite("z_cut", &SVT22MagneticField::z_cut)


#if autodiff_FOUND
        .def_readwrite("active_diff", &SVT22MagneticField::active_diff)
        .def_readonly("all_diff", &SVT22MagneticField::all_diff)

        .def("derivative", [](SVT22MagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif

        .def("at_position", [](SVT22MagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif