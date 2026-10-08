#include "ImagineModels/Han.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_han(py::module_ &m) {
    bind_regular_model<HanMagneticField>(m, "HanMagneticField")
        .def_readwrite("R_min", &HanMagneticField::R_min)
        .def_readwrite("R_max", &HanMagneticField::R_max)
        .def_readwrite("R_s", &HanMagneticField::R_s)
        .def(py::init<const std::string &>(), "model"_a)
        .def("set_model", &HanMagneticField::set_model, "model"_a)
        .def_property_readonly("model", &HanMagneticField::model)
        .def_readonly("available_models", &HanMagneticField::available_models);
}
