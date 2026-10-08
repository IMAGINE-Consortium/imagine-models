#include "ImagineModels/TF17.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_tf17(py::module_ &m) {
    bind_regular_model<TF17MagneticField>(m, "TF17MagneticField")
        .def(py::init<const std::string &, const std::string &>(), "disk_model"_a, "halo_model"_a = "C0")
        .def("set_model", &TF17MagneticField::set_model, "disk_model"_a, "halo_model"_a, doc::set_model)
        .def_property_readonly("disk_model", &TF17MagneticField::disk_model)
        .def_property_readonly("halo_model", &TF17MagneticField::halo_model)
        .def_readonly("available_disk_models", &TF17MagneticField::available_disk_models)
        .def_readonly("available_halo_models", &TF17MagneticField::available_halo_models);
}
