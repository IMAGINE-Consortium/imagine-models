#include "../bindings.h"
#include "ImagineModels/Han.h"
#include "../model_bindings.h"

void bind_han(py::module_ &m)
{
    bind_regular_model<HanMagneticField>(m, "HanMagneticField");
}
