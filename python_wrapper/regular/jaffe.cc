#include "ImagineModels/Jaffe.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_jaffe(py::module_ &m) {
    bind_regular_model<JaffeMagneticField>(m, "JaffeMagneticField")
        .def_readwrite("quadruple", &JaffeMagneticField::quadruple)
        .def_readwrite("bss", &JaffeMagneticField::bss)
        .def_readwrite("ring", &JaffeMagneticField::ring)
        .def_readwrite("bar", &JaffeMagneticField::bar)
        .def_readwrite("arm_num", &JaffeMagneticField::arm_num);
}
