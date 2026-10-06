#ifndef SUNWRAPPER_H
#define SUNWRAPPER_H

#include "ImagineModels/Sun.h"
#include "model_bindings.h"

void Sun2008(py::module_ &m)
{
    bind_regular_model<SunMagneticField>(m, "SunMagneticField");
}

#endif
