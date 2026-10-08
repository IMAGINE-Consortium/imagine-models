#include "ImagineModels/Fauvet.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_fauvet(py::module_ &m) {
    bind_regular_model<FauvetMagneticField>(m, "FauvetMagneticField")
        .def_readwrite("b_r_max", &FauvetMagneticField::b_r_max)
        .def_readwrite("b_r_min", &FauvetMagneticField::b_r_min);
}
