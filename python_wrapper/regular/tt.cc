#include "ImagineModels/TT.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_tt(py::module_ &m) {
    bind_regular_model<TTMagneticField>(m, "TTMagneticField")
        .def_readwrite("b_r_max", &TTMagneticField::b_r_max)
        .def_readwrite("b_r_min", &TTMagneticField::b_r_min);
}
