 #include <pybind11/pybind11.h>

#include "RandomJF12.h"

namespace py = pybind11;
using namespace pybind11::literals;

void RandomJF12(py::module_ &m) {

    py::class_<JF12RandomField, RandomVectorField>(m, "JF12RandomField")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readonly("regular_base", &JF12RandomField::regular_base)

        .def_readwrite("apply_spectrum", &JF12RandomField::apply_spectrum)
        .def_readwrite("clean_divergence", &JF12RandomField::clean_divergence)

        .def_readwrite("spectral_offset", &JF12RandomField::spectral_offset)
        .def_readwrite("spectral_slope", &JF12RandomField::spectral_slope)

        .def_readwrite("b0_1", &JF12RandomField::b0_1)
        .def_readwrite("b0_2", &JF12RandomField::b0_2)
        .def_readwrite("b0_3", &JF12RandomField::b0_3)
        .def_readwrite("b0_4", &JF12RandomField::b0_4)
        .def_readwrite("b0_5", &JF12RandomField::b0_5)
        .def_readwrite("b0_6", &JF12RandomField::b0_6)
        .def_readwrite("b0_7", &JF12RandomField::b0_7)
        .def_readwrite("b0_8", &JF12RandomField::b0_8)
        .def_readwrite("b0_int", &JF12RandomField::b0_int)
        .def_readwrite("b0_halo", &JF12RandomField::b0_halo)
        .def_readwrite("r0_halo", &JF12RandomField::r0_halo)
        .def_readwrite("z0_halo", &JF12RandomField::z0_halo)
        .def_readwrite("z0_spiral", &JF12RandomField::z0_spiral)
        .def_readwrite("rho_GC", &JF12RandomField::rho_GC)
        .def_readwrite("Rmax", &JF12RandomField::Rmax)
        .def_readwrite("anisotropy_rho", &JF12RandomField::anisotropy_rho)

        .def("spatial_profile", &JF12RandomField::spatial_profile, "x"_a, "y"_a, "z"_a, py::return_value_policy::move);
}