#include "ImagineModels/Stanev.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_stanev(py::module_ &m) {
    bind_regular_model<StanevBSSMagneticField>(m, "StanevBSSMagneticField")
        .def_readwrite("b_r_max", &StanevBSSMagneticField::b_r_max)
        .def_readwrite("b_r_min", &StanevBSSMagneticField::b_r_min);
}
