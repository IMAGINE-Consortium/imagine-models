#include "../bindings.h"
#include "ImagineModels/Han.h"
#include "../model_bindings.h"

void bind_han(py::module_ &m)
{
    bind_regular_model<HanMagneticField>(m, "HanMagneticField")
        .def_readwrite("R_min", &HanMagneticField::R_min)
        .def_readwrite("R_max", &HanMagneticField::R_max)
        .def_readwrite("R_s", &HanMagneticField::R_s);
}
