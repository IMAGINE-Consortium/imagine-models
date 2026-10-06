#ifndef JAFFEWRAPPER_H
#define JAFFEWRAPPER_H

#include <pybind11/pybind11.h>

#include "ImagineModels/Jaffe.h"

namespace py = pybind11;
using namespace pybind11::literals;

void Jaffe(py::module_ &m)
{
    py::class_<JaffeMagneticField, RegularVectorField>(m, "JaffeMagneticField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("quadruple", &JaffeMagneticField::quadruple)
        .def_readwrite("bss", &JaffeMagneticField::bss)

        .def_readwrite("disk_amp", &JaffeMagneticField::disk_amp)
        .def_readwrite("disk_z0", &JaffeMagneticField::disk_z0)
        .def_readwrite("halo_amp", &JaffeMagneticField::halo_amp)
        .def_readwrite("halo_z0", &JaffeMagneticField::halo_z0)
        .def_readwrite("r_inner", &JaffeMagneticField::r_inner)
        .def_readwrite("r_scale", &JaffeMagneticField::r_scale)
        .def_readwrite("r_peak", &JaffeMagneticField::r_peak)

        .def_readwrite("ring", &JaffeMagneticField::ring)
        .def_readwrite("bar", &JaffeMagneticField::bar)
        // either ring or bar!
        .def_readwrite("ring_amp", &JaffeMagneticField::ring_amp)
        .def_readwrite("ring_r", &JaffeMagneticField::ring_r)
        .def_readwrite("bar_amp", &JaffeMagneticField::bar_amp)
        .def_readwrite("bar_a", &JaffeMagneticField::bar_a)
        .def_readwrite("bar_b", &JaffeMagneticField::bar_b)
        .def_readwrite("bar_phi0", &JaffeMagneticField::bar_phi0)

        .def_readwrite("arm_num", &JaffeMagneticField::arm_num)
        .def_readwrite("arm_r0", &JaffeMagneticField::arm_r0)
        .def_readwrite("arm_z0", &JaffeMagneticField::arm_z0)
        .def_readwrite("arm_phi1", &JaffeMagneticField::arm_phi1)
        .def_readwrite("arm_phi2", &JaffeMagneticField::arm_phi2)
        .def_readwrite("arm_phi3", &JaffeMagneticField::arm_phi3)
        .def_readwrite("arm_phi4", &JaffeMagneticField::arm_phi4)
        .def_readwrite("arm_amp1", &JaffeMagneticField::arm_amp1)
        .def_readwrite("arm_amp2", &JaffeMagneticField::arm_amp2)
        .def_readwrite("arm_amp3", &JaffeMagneticField::arm_amp3)
        .def_readwrite("arm_amp4", &JaffeMagneticField::arm_amp4)
        .def_readwrite("arm_pitch", &JaffeMagneticField::arm_pitch)

        .def_readwrite("comp_c", &JaffeMagneticField::comp_c)
        .def_readwrite("comp_d", &JaffeMagneticField::comp_d)
        .def_readwrite("comp_r", &JaffeMagneticField::comp_r)
        .def_readwrite("comp_p", &JaffeMagneticField::comp_p)
#if IMAGINE_HAS_AUTODIFF
        .def_readwrite("active_diff", &JaffeMagneticField::active_diff)
        .def_readonly("all_diff", &JaffeMagneticField::all_diff)

        .def("derivative", [](JaffeMagneticField &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif

        .def("at_position", [](JaffeMagneticField &self, double x, double y, double z)
            {
            vector f = self.at_position(x, y, z);
            return std::make_tuple(f[0], f[1], f[2]); 
            },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::take_ownership);
}

#endif
