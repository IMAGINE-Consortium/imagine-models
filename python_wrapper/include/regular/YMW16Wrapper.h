#ifndef YMWWRAPPER_H
#define YMWWRAPPER_H

#include <pybind11/pybind11.h>

#include "YMW.h"

namespace py = pybind11;
using namespace pybind11::literals;

void YMW(py::module_ &m) {
    py::class_<YMW16, RegularScalarField>(m, "YMW16")
        .def(py::init<>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def_readwrite("r_warp", &YMW16::t0_r_warp)
        .def_readwrite("r0", &YMW16::r0)
        .def_readwrite("t0_gamma_w", &YMW16::t0_gamma_w)
        // Thick disc
        .def_readwrite("do_thick_disc", &YMW16::do_thick_disc)
        .def_readwrite("t1_ad", &YMW16::t1_ad)
        .def_readwrite("t1_bd", &YMW16::t1_bd)
        .def_readwrite("t1_n1", &YMW16::t1_n1)
        .def_readwrite("t1_h1", &YMW16::t1_h1)
        // Thin disc
        .def_readwrite("do_thin_disc", &YMW16::do_thin_disc)
        .def_readwrite("t2_a2", &YMW16::t2_a2)
        .def_readwrite("t2_b2", &YMW16::t2_b2)
        .def_readwrite("t2_n2", &YMW16::t2_n2)
        .def_readwrite("t2_k2", &YMW16::t2_k2)
        // Spiral arms
        .def_readwrite("do_spiral_arms", &YMW16::do_spiral_arms)
        .def_readwrite("t3_b2s", &YMW16::t3_b2s)
        .def_readwrite("t3_ka", &YMW16::t3_ka)
        .def_readwrite("t3_aa", &YMW16::t3_aa)
        .def_readwrite("t3_ncn", &YMW16::t3_ncn)
        .def_readwrite("t3_wcn", &YMW16::t3_wcn)
        .def_readwrite("t3_thetacn", &YMW16::t3_thetacn)
        .def_readwrite("t3_nsg", &YMW16::t3_nsg)
        .def_readwrite("t3_wsg", &YMW16::t3_wsg)
        .def_readwrite("t3_thetasg", &YMW16::t3_thetasg)
        .def_readwrite("t3_rmin", &YMW16::t3_rmin)
        .def_readwrite("t3_phimin", &YMW16::t3_phimin)
        .def_readwrite("t3_tpitch", &YMW16::t3_tpitch)
        .def_readwrite("t3_narm", &YMW16::t3_narm)
        .def_readwrite("t3_warm", &YMW16::t3_warm)
        // Galactic center
        .def_readwrite("do_galactic_center", &YMW16::do_galactic_center)
        .def_readwrite("t4_ngc", &YMW16::t4_ngc)
        .def_readwrite("t4_agc", &YMW16::t4_agc)
        .def_readwrite("t4_hgc", &YMW16::t4_hgc)
        // Gum
        .def_readwrite("do_gum", &YMW16::do_gum)
        .def_readwrite("t5_kgn", &YMW16::t5_kgn)
        .def_readwrite("t5_ngn", &YMW16::t5_ngn)
        .def_readwrite("t5_wgn", &YMW16::t5_wgn)
        .def_readwrite("t5_agn", &YMW16::t5_agn)
        // Local bubble
        .def_readwrite("do_local_bubble", &YMW16::do_local_bubble)
        .def_readwrite("t6_j_lb", &YMW16::t6_j_lb)
        .def_readwrite("t6_nlb1", &YMW16::t6_nlb1)
        .def_readwrite("t6_detlb1", &YMW16::t6_detlb1)
        .def_readwrite("t6_wlb1", &YMW16::t6_wlb1)
        .def_readwrite("t6_hlb1", &YMW16::t6_hlb1)
        .def_readwrite("t6_thetalb1", &YMW16::t6_thetalb1)
        .def_readwrite("t6_nlb2", &YMW16::t6_nlb2)
        .def_readwrite("t6_detlb2", &YMW16::t6_detlb2)
        .def_readwrite("t6_wlb2", &YMW16::t6_wlb2)
        .def_readwrite("t6_hlb2", &YMW16::t6_hlb2)
        .def_readwrite("t6_thetalb2", &YMW16::t6_thetalb2)
        // Loop
        .def_readwrite("do_loop", &YMW16::do_loop)
        .def_readwrite("t7_nli", &YMW16::t7_nli)
        .def_readwrite("t7_rli", &YMW16::t7_rli)
        .def_readwrite("t7_wli", &YMW16::t7_wli)
        .def_readwrite("t7_detthetali", &YMW16::t7_detthetali)
        .def_readwrite("t7_thetali", &YMW16::t7_thetali)
#if autodiff_FOUND
        .def_readwrite("active_diff", &YMW16::active_diff)
        .def_readonly("all_diff", &YMW16::all_diff)

        .def("derivative", [](YMW16 &self, double x, double y, double z)
            { return self.derivative(x, y, z); },
            "x"_a, "y"_a, "z"_a, py::return_value_policy::move)
#endif

        .def("at_position", &YMW16::at_position, "x"_a, "y"_a, "z"_a, py::return_value_policy::move);
        


}

#endif