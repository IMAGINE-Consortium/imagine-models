#include "../bindings.h"
#include "ImagineModels/Archimedes.h"
#include "../model_bindings.h"

void bind_archimedes(py::module_ &m)
{
    bind_regular_model<ArchimedeanMagneticField>(m, "ArchimedeanMagneticField");
}
