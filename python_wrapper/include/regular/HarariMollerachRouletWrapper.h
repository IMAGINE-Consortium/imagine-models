#ifndef HARARIMOLLERACHROULETWRAPPER_H
#define HARARIMOLLERACHROULETWRAPPER_H

#include "ImagineModels/HarariMollerachRoulet.h"
#include "model_bindings.h"

void HarariMollerachRoulet(py::module_ &m)
{
    bind_regular_model<HMRMagneticField>(m, "HMRMagneticField")
        .def_readwrite("b_r_max", &HMRMagneticField::b_r_max);
}

#endif
