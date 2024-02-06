#include <cmath>
#include "hamunits.h"
#include "StanevBSS.h"

#include "helpers.h"

// https://arxiv.org/abs/astro-ph/9607086, implementation from Hammurabi (old). Implemented is the bisymmetric model
vector StanevBSSMagneticField::_at_position(const double &x, const double &y, const double &z, const StanevBSSMagneticField &p) const
{

    vector B_vec3{{0, 0, 0}};
    const double r = sqrt(x * x + y * y);
    const double phi = atan2(y, x);

    if (r > b_r_max)
    {
        return B_vec3;
    }

    auto phi_prime = p.b_phi0 - phi; // PHIprime running clock-wise from neg. x-axis
    auto beta = 1. / tan(p.b_p * (M_PI / 180.));

    auto B_0 = 3 * p.b_Rsun / b_r_min;
    if (r > b_r_min)
    {
        B_0 = 3 * p.b_Rsun / r;
    }
    

    auto z_0 = p.b_z01;
    if (std::abs(z) > p.b_z0_border)
    {
        z_0 = p.b_z02;
    }
    // eq. 1, 3, 4
    // minus sign before abs(z) added in eq. 4 -> would make no sense otherwise... 
    vector B_cyl{{B_0 * cos(phi_prime - beta * log(r / p.b_r0)) * sin(p.b_p * (M_PI / 180.)) * exp(-std::abs(z) / z_0),
                  -B_0 * cos(phi_prime - beta * log(r / p.b_r0)) * cos(p.b_p * (M_PI / 180.)) * exp(-std::abs(z) / z_0),
                  0.}};

    B_vec3 = Cyl2Cart<vector>(phi, B_cyl);
    return B_vec3;
}

#if autodiff_FOUND

Eigen::MatrixXd StanevBSSMagneticField::_jac(const double &x, const double &y, const double &z, StanevBSSMagneticField &p) const
{
    vector out;
    Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, StanevBSSMagneticField &_p)
                                          { return _p._at_position(_x, _y, _z, _p); },
                                          ad::wrt(p.b_Rsun, p.b_z01, p.b_z02, p.b_z0_border, p.b_r0, p.b_p, p.b_phi0), ad::at(x, y, z, p), out);
    return _filter_diff(_deriv);
}

#endif