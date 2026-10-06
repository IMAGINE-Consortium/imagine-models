#include "../bindings.h"
#include "ImagineModelsRandom/LogNormal.h"

void bind_lognormal(py::module_ &m)
{
    py::class_<LogNormalScalarField, RandomScalarField>(m, "LogNormalScalarField")
        .def(py::init<>())

        .def_readwrite("apply_spectrum", &LogNormalScalarField::apply_spectrum)

        .def_readwrite("log_mean", &LogNormalScalarField::log_mean)

        .def_readwrite("spectral_offset", &LogNormalScalarField::spectral_offset)
        .def_readwrite("spectral_slope", &LogNormalScalarField::spectral_slope);
}
