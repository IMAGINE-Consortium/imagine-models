#ifndef HELIXWRAPPER_H
#define HELIXWRAPPER_H

#include "ImagineModels/Helix.h"
#include "model_bindings.h"

void Helix(py::module_ &m)
{
    bind_regular_model<HelixMagneticField>(m, "HelixMagneticField")
        .def_readwrite("rmin", &HelixMagneticField::rmin)
        .def_readwrite("rmax", &HelixMagneticField::rmax);
}

#endif
