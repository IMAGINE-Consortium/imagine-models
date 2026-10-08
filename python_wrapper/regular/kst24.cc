#include "ImagineModels/KST24.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_kst24(py::module_ &m) {
    bind_regular_model<KST24MagneticField>(m, "KST24MagneticField")
        .def_readwrite("arm_rmin", &KST24MagneticField::arm_rmin)
        .def_readwrite("arm_rmax", &KST24MagneticField::arm_rmax)
        .def_readwrite("arm_zmax", &KST24MagneticField::arm_zmax)
        .def_readwrite("arm_widening", &KST24MagneticField::arm_widening)
        .def_readwrite("spiral_a", &KST24MagneticField::spiral_a)
        .def_readwrite("width_reference_r", &KST24MagneticField::width_reference_r)
        .def_readwrite("width_max", &KST24MagneticField::width_max)
        .def_readwrite("torus_rmin", &KST24MagneticField::torus_rmin)
        .def_readwrite("X_rmin", &KST24MagneticField::X_rmin)
        .def_readwrite("X_zmax", &KST24MagneticField::X_zmax);
}
