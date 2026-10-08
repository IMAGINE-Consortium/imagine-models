#include "../bindings.h"
#include "ImagineModelsRandom/UF26Random.h"

void bind_uf26_random(py::module_ &m) {
    py::class_<UF26RandomField, RandomVectorField>(m, "UF26RandomField")
        .def(py::init<const std::string &>(), "model"_a = "expDisk")
        .def("set_model", &UF26RandomField::set_model, "model"_a)
        .def_property_readonly("model", &UF26RandomField::model)
        .def_readonly("available_models", &UF26RandomField::available_models)
        .def_readwrite("spectral_offset", &UF26RandomField::spectral_offset)
        .def_readwrite("spectral_slope", &UF26RandomField::spectral_slope)
        .def_readwrite("b_disk", &UF26RandomField::b_disk)
        .def_readwrite("z_disk", &UF26RandomField::z_disk)
        .def_readwrite("l_r", &UF26RandomField::l_r)
        .def_readwrite("r_c", &UF26RandomField::r_c)
        .def_readwrite("b_ring", &UF26RandomField::b_ring)
        .def_readwrite("r_ring", &UF26RandomField::r_ring)
        .def_readwrite("w_ring", &UF26RandomField::w_ring)
        .def_readwrite("z_ring", &UF26RandomField::z_ring)
        .def_readwrite("r_max", &UF26RandomField::r_max)
        .def_readwrite("w_max", &UF26RandomField::w_max)
        .def_readwrite("r_sun", &UF26RandomField::r_sun);
}
