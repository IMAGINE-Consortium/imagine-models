#ifndef HANWRAPPER_H
#define HANWRAPPER_H

#include "ImagineModels/Han.h"
#include "model_bindings.h"

void Han2018(py::module_ &m)
{
    bind_regular_model<HanMagneticField>(m, "HanMagneticField");
}

#endif
