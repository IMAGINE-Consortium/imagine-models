#include "../bindings.h"
#include "ImagineModelsRandom/JF12Random.h"

void bind_jf12_random(py::module_ &m) {
    py::class_<JF12RandomField, RandomVectorField>(m, "JF12RandomField")
        .def(py::init<const std::string &>(), "model"_a = "JF12")
        .def("set_model", &JF12RandomField::set_model, "model"_a, doc::set_model)
        .def_property_readonly("model", &JF12RandomField::model)
        .def_readonly("available_models", &JF12RandomField::available_models)
        .def_readwrite("arm_shift", &JF12RandomField::arm_shift)

        .def_readonly("regular_base", &JF12RandomField::regular_base)

        .def_readwrite("spectral_offset", &JF12RandomField::spectral_offset)
        .def_readwrite("spectral_slope", &JF12RandomField::spectral_slope)
        .def_readwrite("f_iso", &JF12RandomField::f_iso)
        .def_readwrite("f_aniso", &JF12RandomField::f_aniso)
        .def_readwrite("beta", &JF12RandomField::beta)

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
        .def_readwrite("Rmax", &JF12RandomField::Rmax);
}
