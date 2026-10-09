#include "../bindings.h"
#include "ImagineModelsRandom/SunRandom.h"

void bind_sun_random(py::module_ &m) {
    py::class_<SunRandomField, RandomVectorField>(m, "SunRandomField")
        .def(py::init<const std::string &>(), "model"_a = "Sun10")
        .def("set_model", &SunRandomField::set_model, "model"_a, doc::set_model)
        .def_property_readonly("model", &SunRandomField::model)
        .def_readonly("available_models", &SunRandomField::available_models)
        .def_readwrite("spectral_offset", &SunRandomField::spectral_offset)
        .def_readwrite("spectral_slope", &SunRandomField::spectral_slope)
        .def_readwrite("b_iso", &SunRandomField::b_iso)
        .def_readwrite("r_sun", &SunRandomField::r_sun)
        .def_readwrite("r0", &SunRandomField::r0)
        .def_readwrite("h_disk", &SunRandomField::h_disk)
        .def_readwrite("h_halo", &SunRandomField::h_halo)
        .def_readwrite("f_disk", &SunRandomField::f_disk);
}
