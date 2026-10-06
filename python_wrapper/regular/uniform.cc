#include "../bindings.h"
#include "ImagineModels/Uniform.h"
#include "../model_bindings.h"

void bind_uniform(py::module_ &m)
{
    bind_regular_model<UniformMagneticField>(m, "UniformMagneticField");
    bind_regular_model<UniformDensityField>(m, "UniformDensityField");
}
