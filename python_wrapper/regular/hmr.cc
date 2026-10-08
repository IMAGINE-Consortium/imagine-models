#include "../bindings.h"
#include "../model_bindings.h"
#include "ImagineModels/HarariMollerachRoulet.h"

void bind_hmr(py::module_ &m) {
    bind_regular_model<HMRMagneticField>(m, "HMRMagneticField").def_readwrite("b_r_max", &HMRMagneticField::b_r_max);
}
