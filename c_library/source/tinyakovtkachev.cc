#include <cmath>

#include "TinyakovTkachev.h"
#include "helpers.h"

vector TTMagneticField::_at_position(const double &x, const double &y, const double &z, const TTMagneticField &p) const
{

    double r = sqrt(x * x + y * y);
    double phi = atan2(y, x);

    vector B_vec3{{0, 0, 0}};
    if (r > b_r_max)
    {
        return B_vec3;
    }

    auto pitch = p.b_p * (M_PI / 180.);

    auto beta = 1. / tan(pitch);

    auto phase = (beta * log(1. + p.b_d / p.b_Rsun)) - M_PI / 2.;  // eq. 2
    //  double epsilon0=(b5_Rsun+b5_d)*exp(-(M_PI/2.)*tan(b5_p)); <-- hammurabi comment

    double sign = -1.;
    if (z < 0) {  // halo field anti-parallel above/below disk, see eq. 5
        sign = 1.;
    }

    auto f_z = sign * exp(-(std::abs(z) / p.b_z0));  // eq. 5 -> dipole model implememted

    // there is a factor 1/cos(phase) difference between
    // the original TT and Kachelriess. <-- hammurabi comment
    // double b_r=b_b0*(b_Rsun/(r));
    auto b_r = p.b_b0 * (p.b_Rsun / (b_r_min * cos(phase)));
    if (r > b_r_min) {
        b_r = p.b_b0 * (p.b_Rsun / (r * cos(phase)));
    }

    // if(r<b5_r_min) {b_r=b5_b0*(b5_Rsun/(b5_r_min));} <-- hammurabi comment
    else {
        b_r = p.b_b0 * (p.b_Rsun / (b_r_min * cos(phase)));  // eq. 3
    }

    // added + pi/2 because coordinate system / Earth position is rotated by pi/2 in TT paper
    auto B_r_phi = b_r * cos(-phi - beta * log(r / p.b_Rsun) + phase + M_PI / 2.);  // eq. 1
    //  if(r<1.e-26){B_r_phi = b_r*std::cos(phi+phase);} <-- hammurabi comment

    // B-field in cylindrical coordinates: <-- hammurabi comment
    vector B_cyl{{B_r_phi * sin(pitch) * f_z,
                  B_r_phi * cos(pitch) * f_z,
                  0.}};

    B_vec3 = Cyl2Cart<vector>(phi, B_cyl);

    std::cout << "Bx: " << B_vec3[0] << std::endl;
    std::cout << "By: " << B_vec3[1] << std::endl;
    std::cout << "Bz: " << B_vec3[2] << std::endl;

    return B_vec3;
}

#if autodiff_FOUND

Eigen::MatrixXd TTMagneticField::_jac(const double &x, const double &y, const double &z, TTMagneticField &p) const
{
    vector out;
    Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, TTMagneticField &_p)
                                          { return _p._at_position(_x, _y, _z, _p); },
                                          ad::wrt(p.b_Rsun, p.b_b0, p.b_d, p.b_z0, p.b_p), ad::at(x, y, z, p), out);
    return _filter_diff(_deriv);
}

#endif