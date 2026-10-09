#include "ImagineModels/Sun.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_sun(py::module_ &m) {
    bind_regular_model<SunMagneticField>(m, "SunMagneticField")
        .def(py::init<const std::string &>(), "model"_a)
        .def("set_model", &SunMagneticField::set_model, "model"_a, doc::set_model)
        .def_property_readonly("model", &SunMagneticField::model)
        .def_readonly("available_models", &SunMagneticField::available_models);
}
