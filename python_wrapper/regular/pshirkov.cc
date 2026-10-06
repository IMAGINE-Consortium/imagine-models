#include "../bindings.h"
#include "ImagineModels/Pshirkov.h"
#include "../model_bindings.h"

void bind_pshirkov(py::module_ &m)
{
    bind_regular_model<PshirkovMagneticField>(m, "PshirkovMagneticField")
        .def_readwrite("useASS", &PshirkovMagneticField::useASS)
        .def_readwrite("useBSS", &PshirkovMagneticField::useBSS)
        .def_readwrite("useHalo", &PshirkovMagneticField::useHalo)
        .def_readwrite("R_c", &PshirkovMagneticField::R_c);
}
