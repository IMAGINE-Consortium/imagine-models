#ifndef STANEVBSSWRAPPER_H
#define STANEVBSSWRAPPER_H

#include "ImagineModels/StanevBSS.h"
#include "model_bindings.h"

void StanevBSS(py::module_ &m)
{
    bind_regular_model<StanevBSSMagneticField>(m, "StanevBSSMagneticField")
        .def_readwrite("b_r_max", &StanevBSSMagneticField::b_r_max)
        .def_readwrite("b_r_min", &StanevBSSMagneticField::b_r_min);
}

#endif
