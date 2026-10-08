#include "ImagineModels/Sun.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_sun(py::module_ &m) {
    bind_regular_model<SunMagneticField>(m, "SunMagneticField");
}
