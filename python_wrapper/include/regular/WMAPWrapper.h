#ifndef WMAPWRAPPER_H
#define WMAPWRAPPER_H

#include "ImagineModels/WMAP.h"
#include "model_bindings.h"

void WMAP(py::module_ &m)
{
    bind_regular_model<WMAPMagneticField>(m, "WMAPMagneticField")
        .def_readwrite("b_r_max", &WMAPMagneticField::b_r_max)
        .def_readwrite("b_r_min", &WMAPMagneticField::b_r_min)
        .def_readwrite("b_anti", &WMAPMagneticField::anti);
}

#endif
