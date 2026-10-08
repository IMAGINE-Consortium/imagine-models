#include "ImagineModels/Uniform.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_uniform(py::module_ &m) {
    bind_regular_model<UniformMagneticField>(m, "UniformMagneticField");
    bind_regular_model<UniformDensityField>(m, "UniformDensityField");
}
