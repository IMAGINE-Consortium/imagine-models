#include "ImagineModels/JF12.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_jf12(py::module_ &m) {
    bind_regular_model<JF12MagneticField>(m, "JF12MagneticField")
        .def_readwrite("do_halo", &JF12MagneticField::do_halo)
        .def_readwrite("do_X", &JF12MagneticField::do_X)
        .def(py::init<const std::string &>(), "model"_a)
        .def("set_model", &JF12MagneticField::set_model, "model"_a)
        .def_property_readonly("model", &JF12MagneticField::model)
        .def_readonly("available_models", &JF12MagneticField::available_models)
        .def_readwrite("arm_shift", &JF12MagneticField::arm_shift);
}
