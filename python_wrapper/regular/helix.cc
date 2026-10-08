#include "ImagineModels/Helix.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_helix(py::module_ &m) {
    bind_regular_model<HelixMagneticField>(m, "HelixMagneticField")
        .def_readwrite("rmin", &HelixMagneticField::rmin)
        .def_readwrite("rmax", &HelixMagneticField::rmax);
}
