#include "ImagineModels/SVT22.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_svt22(py::module_ &m) {
    bind_regular_model<SVT22MagneticField>(m, "SVT22").def_readwrite("do_halo", &SVT22MagneticField::do_halo);
}
