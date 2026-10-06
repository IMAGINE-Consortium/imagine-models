#ifndef SVT22WRAPPER_H
#define SVT22WRAPPER_H

#include "ImagineModels/SVT22.h"
#include "model_bindings.h"

void SVT22(py::module_ &m)
{
    bind_regular_model<SVT22MagneticField>(m, "SVT22");
}

#endif
