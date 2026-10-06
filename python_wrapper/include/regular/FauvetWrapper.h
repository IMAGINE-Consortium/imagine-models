#ifndef FAUVETWRAPPER_H
#define FAUVETWRAPPER_H

#include "ImagineModels/Fauvet.h"
#include "model_bindings.h"

void Fauvet(py::module_ &m)
{
    bind_regular_model<FauvetMagneticField>(m, "FauvetMagneticField")
        .def_readwrite("b_r_max", &FauvetMagneticField::b_r_max)
        .def_readwrite("b_r_min", &FauvetMagneticField::b_r_min);
}

#endif
