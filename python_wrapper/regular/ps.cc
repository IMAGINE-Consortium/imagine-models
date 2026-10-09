#include "ImagineModels/PS.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_ps(py::module_ &m) {
    bind_regular_model<PSMagneticField>(m, "PSMagneticField")
        .def_readwrite("do_disk", &PSMagneticField::do_disk)
        .def_readwrite("do_halo", &PSMagneticField::do_halo)
        .def_readwrite("do_dipole", &PSMagneticField::do_dipole)
        .def_readwrite("b_r_max", &PSMagneticField::b_r_max)
        .def_readwrite("b_r_min", &PSMagneticField::b_r_min)
        .def_readwrite("h_R0", &PSMagneticField::h_R0)
        .def_readwrite("d_r_core", &PSMagneticField::d_r_core)
        .def_readwrite("d_b_core", &PSMagneticField::d_b_core);
}
