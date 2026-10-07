#include "../bindings.h"
#include "ImagineModels/Pshirkov.h"
#include "../model_bindings.h"

void bind_pshirkov(py::module_ &m)
{
    bind_regular_model<PshirkovMagneticField>(m, "PshirkovMagneticField")
        .def(py::init<const std::string &>(), "model"_a)
        .def("set_model", &PshirkovMagneticField::set_model, "model"_a)
        .def_property_readonly("model", &PshirkovMagneticField::model)
        .def_readonly("available_models", &PshirkovMagneticField::available_models)
        .def_readwrite("useDisk", &PshirkovMagneticField::useDisk)
        .def_readwrite("useHalo", &PshirkovMagneticField::useHalo)
        .def_readwrite("R_c", &PshirkovMagneticField::R_c);
}
