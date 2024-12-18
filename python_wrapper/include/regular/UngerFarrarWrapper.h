#ifndef UF24WRAPPER_H
#define UF24WRAPPER_H

#include <pybind11/pybind11.h>

#include "UngerFarrar.h"

void UF24(py::module_ &m)
{
    py::class_<UFMagneticField, RegularVectorField>(m, "UFMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("fDiskB1", &UFMagneticField::fDiskB1)
        .def_readwrite("fDiskB2", &UFMagneticField::fDiskB2)
        .def_readwrite("fDiskB3", &UFMagneticField::fDiskB3)
        .def_readwrite("fDiskH", &UFMagneticField::fDiskH)
        .def_readwrite("fDiskPhase1", &UFMagneticField::fDiskPhase1)
        .def_readwrite("fDiskPhase2", &UFMagneticField::fDiskPhase2)
        .def_readwrite("fDiskPhase3", &UFMagneticField::fDiskPhase3)
        .def_readwrite("fDiskPitch", &UFMagneticField::fDiskPitch)
        .def_readwrite("fDiskW", &UFMagneticField::fDiskW)
        .def_readwrite("fPoloidalA", &UFMagneticField::fPoloidalA)
        .def_readwrite("fPoloidalB", &UFMagneticField::fPoloidalB)
        .def_readwrite("fPoloidalP", &UFMagneticField::fPoloidalP)
        .def_readwrite("fPoloidalR", &UFMagneticField::fPoloidalR)
        .def_readwrite("fPoloidalW", &UFMagneticField::fPoloidalW)
        .def_readwrite("fPoloidalZ", &UFMagneticField::fPoloidalZ)
        .def_readwrite("fSpurCenter", &UFMagneticField::fSpurCenter)
        .def_readwrite("fSpurLength", &UFMagneticField::fSpurLength)
        .def_readwrite("fSpurWidth", &UFMagneticField::fSpurWidth)
        .def_readwrite("fStriation", &UFMagneticField::fStriation)
        .def_readwrite("fToroidalBN", &UFMagneticField::fToroidalBN)
        .def_readwrite("fToroidalBS", &UFMagneticField::fToroidalBS)
        .def_readwrite("fToroidalR", &UFMagneticField::fToroidalR)
        .def_readwrite("fToroidalW", &UFMagneticField::fToroidalW)
        .def_readwrite("fToroidalZ", &UFMagneticField::fToroidalZ)
        .def_readwrite("fTwistingTime", &UFMagneticField::fTwistingTime)

        .def_readwrite("activeModel", &UFMagneticField::activeModel)
        .def_readonly("possibleModels", &UFMagneticField::possibleModels)
        .def_readonly("all_parameters", &UFMagneticField::all_parameters)

#if autodiff_FOUND
        .def_readwrite("active_diff", &UFMagneticField::active_diff)
        .def_readonly("all_diff", &UFMagneticField::all_diff)

        .def(
            "derivative", [](UFMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif
        .def(
            "at_position", [](UFMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership)

        .def("set_parameters", [](UFMagneticField &self, std::string model_type) {
            self.set_parameters(model_type); 
        });
}

#endif