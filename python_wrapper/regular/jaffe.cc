#include "ImagineModels/Jaffe.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_jaffe(py::module_ &m) {
    bind_regular_model<JaffeMagneticField>(m, "JaffeMagneticField")
        .def(py::init<const std::string &>(), "model"_a)
        .def("set_model", &JaffeMagneticField::set_model, "model"_a, doc::set_model)
        .def_property_readonly("model", &JaffeMagneticField::model)
        .def_readonly("available_models", &JaffeMagneticField::available_models)
        .def_readwrite("quadruple", &JaffeMagneticField::quadruple)
        .def_readwrite("bss", &JaffeMagneticField::bss)
        .def_readwrite("ring", &JaffeMagneticField::ring)
        .def_readwrite("bar", &JaffeMagneticField::bar)
        .def_readwrite("arm_num", &JaffeMagneticField::arm_num)
        .def_readwrite("hammurabi_v3", &JaffeMagneticField::hammurabi_v3)
        .def_readwrite("r_max", &JaffeMagneticField::r_max);
}
