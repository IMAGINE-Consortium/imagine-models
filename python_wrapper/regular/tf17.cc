#include "../bindings.h"
#include "ImagineModels/TF17.h"
#include "../model_bindings.h"

void bind_tf17(py::module_ &m)
{
    bind_regular_model<TFMagneticField>(m, "TFMagneticField")
        .def(py::init<const std::string &, const std::string &>(), "disk_model"_a, "halo_model"_a = "C0")
        .def("set_model", &TFMagneticField::set_model, "disk_model"_a, "halo_model"_a)
        .def_property_readonly("disk_model", &TFMagneticField::disk_model)
        .def_property_readonly("halo_model", &TFMagneticField::halo_model)
        .def_readonly("available_disk_models", &TFMagneticField::available_disk_models)
        .def_readonly("available_halo_models", &TFMagneticField::available_halo_models);
}
