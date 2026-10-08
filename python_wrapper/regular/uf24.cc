#include "ImagineModels/UF24.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_uf24(py::module_ &m) {
    bind_regular_model<UFMagneticField>(m, "UFMagneticField")
        .def(py::init<const std::string &>(), "model"_a)
        .def("set_model", &UFMagneticField::set_model, "model"_a)
        .def_property_readonly("model", &UFMagneticField::model)
        .def_readwrite("fMaxRadius", &UFMagneticField::fMaxRadius)
        .def_readonly("available_models", &UFMagneticField::available_models)
        .def_readonly("all_parameters", &UFMagneticField::all_parameters);
}
