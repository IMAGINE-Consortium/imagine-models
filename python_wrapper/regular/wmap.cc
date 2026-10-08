#include "ImagineModels/WMAP.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_wmap(py::module_ &m) {
    bind_regular_model<WMAPMagneticField>(m, "WMAPMagneticField")
        .def_readwrite("b_r_max", &WMAPMagneticField::b_r_max)
        .def_readwrite("b_r_min", &WMAPMagneticField::b_r_min)
        .def_readwrite("b_anti", &WMAPMagneticField::anti);
}
