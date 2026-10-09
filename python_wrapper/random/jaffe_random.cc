#include "../bindings.h"
#include "ImagineModelsRandom/JaffeRandom.h"

void bind_jaffe_random(py::module_ &m) {
    py::class_<JaffeRandomField, RandomVectorField>(m, "JaffeRandomField")
        .def(py::init<const std::string &>(), "model"_a = "Jaffe13")
        .def("set_model", &JaffeRandomField::set_model, "model"_a, doc::set_model)
        .def_property_readonly("model", &JaffeRandomField::model)
        .def_readonly("available_models", &JaffeRandomField::available_models)
        .def_readonly("regular_base", &JaffeRandomField::regular_base)
        .def_readwrite("spectral_offset", &JaffeRandomField::spectral_offset)
        .def_readwrite("spectral_slope", &JaffeRandomField::spectral_slope)
        .def_readwrite("b_rms", &JaffeRandomField::b_rms)
        .def_readwrite("h_rms", &JaffeRandomField::h_rms)
        .def_readwrite("r_grf", &JaffeRandomField::r_grf)
        .def_readwrite("f_ord", &JaffeRandomField::f_ord);
}
