#include "../bindings.h"
#include "ImagineModelsRandom/Orlando26Random.h"

void bind_orlando26_random(py::module_ &m) {
    py::class_<Orlando26RandomField, RandomVectorField>(m, "Orlando26RandomField")
        .def(py::init<const std::string &>(), "model"_a = "halo4kpc")
        .def("set_model", &Orlando26RandomField::set_model, "model"_a, doc::set_model)
        .def_property_readonly("model", &Orlando26RandomField::model)
        .def_readonly("available_models", &Orlando26RandomField::available_models)
        .def_readonly("regular_base", &Orlando26RandomField::regular_base)
        .def_readwrite("spectral_offset", &Orlando26RandomField::spectral_offset)
        .def_readwrite("spectral_slope", &Orlando26RandomField::spectral_slope)
        .def_readwrite("b_ran", &Orlando26RandomField::b_ran)
        .def_readwrite("r0_ran", &Orlando26RandomField::r0_ran)
        .def_readwrite("z0_ran", &Orlando26RandomField::z0_ran)
        .def_readwrite("r_sun", &Orlando26RandomField::r_sun)
        .def_readwrite("b_ordered", &Orlando26RandomField::b_ordered);
}
