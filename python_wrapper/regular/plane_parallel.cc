#include "../bindings.h"
#include "../model_bindings.h"
#include "ImagineModels/PlaneParallel.h"

void bind_plane_parallel(py::module_ &m) {
    bind_regular_model<PlaneParallelDensity>(m, "PlaneParallelDensity")
        .def(py::init<const std::string &>(), "model"_a)
        .def("set_model", &PlaneParallelDensity::set_model, "model"_a, doc::set_model)
        .def_property_readonly("model", &PlaneParallelDensity::model)
        .def_readonly("available_models", &PlaneParallelDensity::available_models);
}
