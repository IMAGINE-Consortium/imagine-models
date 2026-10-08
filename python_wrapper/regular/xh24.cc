#include "ImagineModels/XH24.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_xh24(py::module_ &m) {
    bind_regular_model<XH24MagneticField>(m, "XH24MagneticField").def_readwrite("r_max", &XH24MagneticField::r_max);
}
