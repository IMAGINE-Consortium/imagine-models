#include <pybind11/pybind11.h>

#include "EnsslinSteininger.h"

namespace py = pybind11;
using namespace pybind11::literals;

void EnsslinSteininger(py::module_ &m) {
    py::class_<ESRandomField, RandomVectorField>(m, "ESRandomField")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def_readwrite("apply_spectrum", &ESRandomField::apply_spectrum)

        .def_readwrite("spectral_offset", &ESRandomField::spectral_offset)
        .def_readwrite("spectral_slope", &ESRandomField::spectral_slope)

        .def_readwrite("r0", &ESRandomField::r0)
        .def_readwrite("z0", &ESRandomField::z0);

}