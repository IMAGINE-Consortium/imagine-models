#include "../bindings.h"
#include "ImagineModelsRandom/GaussianScalar.h"

void bind_gaussian_scalar(py::module_ &m) {
    py::class_<GaussianScalarField, RandomScalarField>(m, "GaussianScalarField")
        .def(py::init<>())
        .def_readwrite("mu", &GaussianScalarField::mu)
        .def_readwrite("sigma", &GaussianScalarField::sigma)
        .def_readwrite("spectral_offset", &GaussianScalarField::spectral_offset)
        .def_readwrite("spectral_slope", &GaussianScalarField::spectral_slope);
}
