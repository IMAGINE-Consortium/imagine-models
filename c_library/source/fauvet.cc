#include <cmath>
#include "hamunits.h"
#include "Fauvet.h"
#include "helpers.h"

// ??????, implementation from Hammurabi (old)

vector FauvetMagneticField::_at_position(const double &x, const double &y, const double &z,  const FauvetMagneticField &p) const
{
    vector B_vec3{{0, 0, 0}};
    const double r = sqrt(x * x + y * y);

    if (r > b_r_max || r < b_r_min)
    {
        return B_vec3;
    }

    double phi = atan2(y, x);
    auto chi_z = p.b_chi0 * (M_PI / 180.) * tanh(z / p.b_z0);
    auto beta = 1. / tan(p.b_p * (M_PI / 180.));

    // B-field in cylindrical coordinates:
    vector B_cyl{{p.b_b0 * cos(phi + beta * log(r / p.b_r0)) * sin(p.b_p * (M_PI / 180.)) * cos(chi_z),
                  -p.b_b0 * cos(phi + beta * log(r / p.b_r0)) * cos(p.b_p * (M_PI / 180.)) * cos(chi_z),
                  p.b_b0 * sin(chi_z)}};

    // Taking into account the halo field
    number h_z1;
    if (std::abs(z) < p.h_z0)
    {
        h_z1 = p.h_z1a;
    }
    else
    {
        h_z1 = p.h_z1b;
    }

    auto hf_piece1 = (h_z1 * h_z1) / (h_z1 * h_z1 + (std::abs(z) - p.h_z0) * (std::abs(z) - p.h_z0));
    auto hf_piece2 = exp(-(r - p.h_r0) / (p.h_r0));

    auto halo_field = p.h_b0 * hf_piece1 * (r / p.b_r0) * hf_piece2;
    B_cyl[1] += halo_field;

    B_vec3 = Cyl2Cart<vector>(phi, B_cyl);
    return B_vec3;
}

#if autodiff_FOUND

Eigen::MatrixXd FauvetMagneticField::_jac(const double &x, const double &y, const double &z, FauvetMagneticField &p) const
{
  vector out;
  Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, FauvetMagneticField &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(p.b_b0, p.b_z0, p.b_r0, p.b_p, p.b_chi0, p.h_b0, p.h_z0, p.h_r0, p.h_z1a, p.h_z1b), ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif

