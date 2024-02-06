#ifndef JAFFE_H
#define JAFFE_H

#include <cmath>
#include <vector>
#include <iostream>

#include "Field.h"
#include "RegularField.h"

class JaffeMagneticField : public RegularVectorField
{
protected:
    vector _at_position(const double &x, const double &y, const double &z, const JaffeMagneticField &p) const;

#if autodiff_FOUND
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, JaffeMagneticField &p) const;
#endif

public:
    using RegularVectorField ::RegularVectorField;

    bool quadruple = false; // quadruple pattern in halo
    bool bss = false;       // bi-symmetric

    number disk_amp = 0.167; // disk amplitude, microG
    number disk_z0 = 0.1;    // disk height scale, kpc
    number halo_amp = 1.38;  // halo amplitude, microG
    number halo_z0 = 3.0;    // halo height scale, kpc
    number r_inner = 0.5;    // inner R scale, kpc
    number r_scale = 20.;    // R scale, kpc
    number r_peak = 0.;      // R peak, kpc

    bool ring = false;       // molecular ring
    bool bar = true;         // elliptical bar
                             // either ring or bar!
    number ring_amp = 0.023; // ring field amplitude, microG
    number ring_r = 5.0;     // ring radius, kpc
    number bar_amp = 0.023;  // bar field amplitude, microG
    number bar_a = 5.0;      // major scale, kpc
    number bar_b = 3.0;      // minor scale, kpc
    number bar_phi0 = 45.0;  // bar major direction

    int arm_num = 4;      // # of spiral arms
    number arm_r0 = 7.1;  // arm ref radius, kpc
    number arm_z0 = 0.1;  // arm heigth scale, kpc
    number arm_phi1 = 70; // arm ref angles, deg
    number arm_phi2 = 160;
    number arm_phi3 = 250;
    number arm_phi4 = 340;
    number arm_amp1 = 2; //  arm field amplitudes, microG
    number arm_amp2 = 0.133;
    number arm_amp3 = -3.78;
    number arm_amp4 = 0.32;
    number arm_pitch = 11.5; // pitch angle, deg

    number comp_c = 0.5; // compress factor
    number comp_d = 0.3; // arm cross-sec scale, kpc
    number comp_r = 12;  // radial cutoff scale, kpc
    number comp_p = 3;   // cutoff power

#if autodiff_FOUND
    const std::set<std::string> all_diff{"disk_amp", "disk_z0", "halo_amp", "halo_z0", "r_inner", "r_scale", "r_peak",
                                         "ring_amp", "ring_r", "bar_amp", "bar_a", "bar_b", "bar_phi0",
                                         "arm_r0", "arm_z0", "arm_phi1", "arm_phi2", "arm_phi3", "arm_phi4",
                                         "arm_amp1", "arm_amp2", "arm_amp3", "arm_amp4", "arm_pitch",
                                         "comp_c", "comp_d", "comp_r", "comp_p"};
    std::set<std::string> active_diff{"disk_amp", "disk_z0", "halo_amp", "halo_z0", "r_inner", "r_scale", "r_peak",
                                      "ring_amp", "ring_r", "bar_amp", "bar_a", "bar_b", "bar_phi0",
                                      "arm_r0", "arm_z0", "arm_phi1", "arm_phi2", "arm_phi3", "arm_phi4",
                                      "arm_amp1", "arm_amp2", "arm_amp3", "arm_amp4", "arm_pitch",
                                      "comp_c", "comp_d", "comp_r", "comp_p"};

    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
    {
        return _jac(x, y, z, *this);
    }
#endif

    vector at_position(const double &x, const double &y, const double &z) const
    {
        return _at_position(x, y, z, *this);
    }

    vector orientation(const double &x, const double &y, const double &z, const JaffeMagneticField &p) const;

    number radial_scaling(const double &x, const double &y, const JaffeMagneticField &p) const;

    std::vector<number> arm_compress(const double &x, const double &y, const double &z, const JaffeMagneticField &p) const;

    std::vector<number> arm_compress_dust(const double &x, const double &y, const double &z, const JaffeMagneticField &p) const;

    std::vector<number> dist2arm(const double &x, const double &y, const JaffeMagneticField &p) const;

    number arm_scaling(const double &z, const JaffeMagneticField &p) const;

    number disk_scaling(const double &z, const JaffeMagneticField &p) const;

    number halo_scaling(const double &z, const JaffeMagneticField &p) const;
};

#endif