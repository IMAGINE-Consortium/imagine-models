#ifndef ARCHIMEDEANWRAPPER_H
#define ARCHIMEDEANWRAPPER_H

#include "ImagineModels/Archimedes.h"
#include "model_bindings.h"

void Archimedes(py::module_ &m)
{
    bind_regular_model<ArchimedeanMagneticField>(m, "ArchimedeanMagneticField");
}

#endif
