#include "../bindings.h"
#include "ImagineModels/YMW.h"
#include "../model_bindings.h"

void bind_ymw16(py::module_ &m)
{
    bind_regular_model<YMW16>(m, "YMW16")
        .def_readwrite("r_warp", &YMW16::t0_r_warp)
        .def_readwrite("t0_gamma_w", &YMW16::t0_gamma_w)
        .def_readwrite("do_thick_disc", &YMW16::do_thick_disc)
        .def_readwrite("do_thin_disc", &YMW16::do_thin_disc)
        .def_readwrite("do_spiral_arms", &YMW16::do_spiral_arms)
        .def_readwrite("t3_rmin", &YMW16::t3_rmin)
        .def_readwrite("t3_phimin", &YMW16::t3_phimin)
        .def_readwrite("t3_tpitch", &YMW16::t3_tpitch)
        .def_readwrite("t3_narm", &YMW16::t3_narm)
        .def_readwrite("t3_warm", &YMW16::t3_warm)
        .def_readwrite("do_galactic_center", &YMW16::do_galactic_center)
        .def_readwrite("do_gum", &YMW16::do_gum)
        .def_readwrite("do_local_bubble", &YMW16::do_local_bubble)
        .def_readwrite("do_loop", &YMW16::do_loop);
}
