#ifndef TF17WRAPPER_H
#define TF17WRAPPER_H

#include <pybind11/pybind11.h>

#include "TF17.h"

void TF17(py::module_ &m)
{
    py::class_<TFMagneticField, RegularVectorField>(m, "TFMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("a_disk", &TFMagneticField::a_disk)
        .def_readwrite("z1_disk", &TFMagneticField::z1_disk)
        .def_readwrite("r1_disk", &TFMagneticField::r1_disk)
        .def_readwrite("B1_disk", &TFMagneticField::B1_disk)
        .def_readwrite("L_disk", &TFMagneticField::L_disk)
        .def_readwrite("phi_star_disk", &TFMagneticField::phi_star_disk)
        .def_readwrite("H_disk", &TFMagneticField::H_disk)
        .def_readwrite("a_halo", &TFMagneticField::a_halo)
        .def_readwrite("z1_halo", &TFMagneticField::z1_halo)
        .def_readwrite("B0_Hs", &TFMagneticField::z1_halo)
        .def_readwrite("B1_halo", &TFMagneticField::B1_halo)
        .def_readwrite("L_halo", &TFMagneticField::L_halo)
        .def_readwrite("phi_star_halo", &TFMagneticField::phi_star_halo)
        .def_readwrite("p_0", &TFMagneticField::p_0)
        .def_readwrite("H_p", &TFMagneticField::H_p)
        .def_readwrite("L_p", &TFMagneticField::L_p)

        .def_readwrite("activeDiskModel", &TFMagneticField::activeDiskModel)
        .def_readwrite("activeHaloModel", &TFMagneticField::activeHaloModel)
        .def_readonly("possibleDiskModels", &TFMagneticField::possibleDiskModels)
        .def_readonly("possibleHaloModels", &TFMagneticField::possibleHaloModels)

#if autodiff_FOUND
        .def_readwrite("active_diff", &TFMagneticField::active_diff)
        .def_readonly("all_diff", &TFMagneticField::all_diff)

        .def(
            "derivative", [](TFMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif
        .def(
            "at_position", [](TFMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership)

        .def("set_params", [](TFMagneticField &self, std::string dtype, std::string htype) {
            self.set_params(dtype, htype); 
        });
}

#endif