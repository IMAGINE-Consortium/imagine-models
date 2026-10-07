#include "../bindings.h"
#include "ImagineModelsRandom/EnsslinSteininger.h"

void bind_es_random(py::module_ &m)
{
    py::class_<ESRandomField, RandomVectorField>(m, "ESRandomField")
        .def(py::init<>())
        .def_readwrite("b0", &ESRandomField::b0)
        .def_readwrite("r0", &ESRandomField::r0)
        .def_readwrite("z0", &ESRandomField::z0)
        .def_readwrite("observer", &ESRandomField::observer)
        .def_readwrite("spectral_offset", &ESRandomField::spectral_offset)
        .def_readwrite("spectral_slope", &ESRandomField::spectral_slope);
}
