#include "../bindings.h"
#include "../model_bindings.h"
#include "ImagineModels/YMW.h"

void bind_ymw16(py::module_ &m) {
    bind_regular_model<YMW16>(m, "YMW16")
        .def_readwrite("max_radius", &YMW16::max_radius)
        .def_readwrite("r_warp", &YMW16::t0_r_warp)
        .def_readwrite("t0_gamma_w", &YMW16::t0_gamma_w)
        .def_readwrite("t0_theta0", &YMW16::t0_theta0)
        .def_readwrite("h0", &YMW16::h0)
        .def_readwrite("h1", &YMW16::h1)
        .def_readwrite("h2", &YMW16::h2)
        .def_readwrite("localbubble_boundary", &YMW16::localbubble_boundary)
        .def_readwrite("Xgc", &YMW16::Xgc)
        .def_readwrite("Ygc", &YMW16::Ygc)
        .def_readwrite("Zgc", &YMW16::Zgc)
        .def_readwrite("t5_lc", &YMW16::t5_lc)
        .def_readwrite("t6_zyl1", &YMW16::t6_zyl1)
        .def_readwrite("t6_zyl2", &YMW16::t6_zyl2)
        .def_readwrite("do_thick_disc", &YMW16::do_thick_disc)
        .def_readwrite("do_thin_disc", &YMW16::do_thin_disc)
        .def_readwrite("do_spiral_arms", &YMW16::do_spiral_arms)
        .def_readwrite("t3_rmin", &YMW16::t3_rmin)
        .def_readwrite("t3_thmin", &YMW16::t3_thmin)
        .def_readwrite("t3_tan_pitch", &YMW16::t3_tan_pitch)
        .def_readwrite("t3_cos_pitch", &YMW16::t3_cos_pitch)
        .def_readwrite("t3_narm", &YMW16::t3_narm)
        .def_readwrite("t3_warm", &YMW16::t3_warm)
        .def_readwrite("do_galactic_center", &YMW16::do_galactic_center)
        .def_readwrite("do_gum", &YMW16::do_gum)
        .def_readwrite("do_local_bubble", &YMW16::do_local_bubble)
        .def_readwrite("do_loop", &YMW16::do_loop);
}
