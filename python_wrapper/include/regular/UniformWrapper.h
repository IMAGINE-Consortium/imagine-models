#ifndef UNIFORMWRAPPER_H
#define UNIFORMWRAPPER_H

#include "ImagineModels/Uniform.h"
#include "model_bindings.h"

void Uniform(py::module_ &m)
{
    bind_regular_model<UniformMagneticField>(m, "UniformMagneticField");
    bind_regular_model<UniformDensityField>(m, "UniformDensityField");
}

#endif
