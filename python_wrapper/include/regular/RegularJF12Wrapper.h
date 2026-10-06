#ifndef REGULARJF12WRAPPER_H
#define REGULARJF12WRAPPER_H

#include "ImagineModels/RegularJF12.h"
#include "model_bindings.h"

void RegularJF12(py::module_ &m)
{
    bind_regular_model<JF12MagneticField>(m, "JF12RegularField")
        .def_readwrite("do_halo", &JF12MagneticField::do_halo)
        .def_readwrite("do_X", &JF12MagneticField::do_X);
}

#endif
