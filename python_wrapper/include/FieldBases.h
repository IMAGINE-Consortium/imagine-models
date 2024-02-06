#ifndef FIELDBASES_H
#define FIELDBASES_H

#include <pybind11/pybind11.h>

#include "regular_trampoline.h"

#if FFTW_FOUND
    #include "random_trampoline.h"
#endif


namespace py = pybind11;
using namespace pybind11::literals;

void FieldBases(py::module_ &m) {
    
    py::class_<Field<vector, std::array<double*, 3>>,  PyVectorFieldBase>(m, "VectorFieldBase");

    py::class_<Field<number, double*>,  PyScalarFieldBase>(m, "ScalarFieldBase");

    #if FFTW_FOUND
        py::class_<RandomField<vector, std::array<double*, 3>>,  PyVectorRandomFieldBase>(m, "VectorRandomFieldBase");

        py::class_<RandomField<number, double*>,  PyScalarRandomFieldBase>(m, "ScalarRandomFieldBase");
    #endif

}

#endif