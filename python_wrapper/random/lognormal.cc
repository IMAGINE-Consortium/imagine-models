#include "ImagineModelsRandom/LogNormal.h"
#include "../bindings.h"

void bind_lognormal(py::module_ &m) {
    py::class_<LogNormalScalarField, RandomScalarField>(m, "LogNormalScalarField")
        .def(py::init<>())
        .def_readwrite("log_mu", &LogNormalScalarField::log_mu)
        .def_readwrite("log_sigma", &LogNormalScalarField::log_sigma)
        .def_readwrite("spectral_offset", &LogNormalScalarField::spectral_offset)
        .def_readwrite("spectral_slope", &LogNormalScalarField::spectral_slope);
}
