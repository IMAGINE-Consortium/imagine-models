#include "ImagineModels/NE2025.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_ne2025(py::module_ &m) {
    bind_regular_model<NE2025>(m, "NE2025")
        .def(py::init<const std::string &>(), "model"_a)
        .def("set_model", &NE2025::set_model, "model"_a, doc::set_model)
        .def_property_readonly("model", &NE2025::model)
        .def_readonly("available_models", &NE2025::available_models);
}
