#ifndef TINYAKOVTKACHEVWRAPPER_H
#define TINYAKOVTKACHEVWRAPPER_H

#include "ImagineModels/TinyakovTkachev.h"
#include "model_bindings.h"

void TinyakovTkachev(py::module_ &m)
{
    bind_regular_model<TTMagneticField>(m, "TTMagneticField")
        .def_readwrite("b_r_max", &TTMagneticField::b_r_max)
        .def_readwrite("b_r_min", &TTMagneticField::b_r_min);
}

#endif
