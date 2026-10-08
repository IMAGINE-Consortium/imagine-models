#include "ImagineModels/UF24.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_uf24(py::module_ &m) {
    bind_regular_model<UF24MagneticField>(m, "UF24MagneticField")
        .def(py::init<const std::string &>(), "model"_a)
        .def("set_model", &UF24MagneticField::set_model, "model"_a, doc::set_model)
        .def_property_readonly("model", &UF24MagneticField::model)
        .def_readwrite("fMaxRadius", &UF24MagneticField::fMaxRadius)
        .def_readonly("available_models", &UF24MagneticField::available_models)
        .def_readonly("all_parameters", &UF24MagneticField::all_parameters);
}
