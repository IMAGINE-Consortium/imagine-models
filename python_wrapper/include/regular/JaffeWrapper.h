#ifndef JAFFEWRAPPER_H
#define JAFFEWRAPPER_H

#include "ImagineModels/Jaffe.h"
#include "model_bindings.h"

void Jaffe(py::module_ &m)
{
    bind_regular_model<JaffeMagneticField>(m, "JaffeMagneticField")
        .def_readwrite("quadruple", &JaffeMagneticField::quadruple)
        .def_readwrite("bss", &JaffeMagneticField::bss)
        .def_readwrite("ring", &JaffeMagneticField::ring)
        .def_readwrite("bar", &JaffeMagneticField::bar)
        .def_readwrite("arm_num", &JaffeMagneticField::arm_num);
}

#endif
