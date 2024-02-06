#include <cmath>
#include "hamunits.h"
#include "WMAP.h"

#include "helpers.h"


// https://iopscience.iop.org/article/10.1086/513699, implementation from Hammurabi (old)

vector WMAPMagneticField::_at_position(const double &x, const double &y, const double &z, const WMAPMagneticField &p) const { 
    
    vector B_vec3{{0, 0, 0}};
    double r = sqrt(x*x + y*y);

    if (r > b_r_max || r < b_r_min) { 
        return B_vec3;
    }

	double phi = atan2(y, x);

    auto psi_r = p.b_psi0*(M_PI/180.) + p.b_psi1*(M_PI/180.) * log(r/p.b_r0);
    auto xsi_z = p.b_xsi0*(M_PI/180.) * tanh(z/p.b_z0);
    
    vector B_cyl{{p.b_b0 * sin(psi_r) * cos(xsi_z),   // eq. 9
                  p.b_b0 * cos(psi_r) * cos(xsi_z), 
                  p.b_b0 * sin(xsi_z) }
                };

    B_vec3 = Cyl2Cart<vector>(phi, B_cyl);
    
    // Antisymmetric, swap the signs.  The way my pitch angle is defined,
    // it seems this has to be swapped this way.  <------ hammurabi comment

    // if (anti && z > 0) {
    //     B_vec3[0] *= (-1.); 
    //     B_vec3[1] *= (-1.); 
    //     B_vec3[2] *= (-1.);
    //     }
    return B_vec3;
}

#if autodiff_FOUND

Eigen::MatrixXd WMAPMagneticField::_jac(const double &x, const double &y, const double &z, WMAPMagneticField &p) const
{
  vector out;
  Eigen::MatrixXd _deriv = ad::jacobian([&](double _x, double _y, double _z, WMAPMagneticField &_p)
                                        { return _p._at_position(_x, _y, _z, _p); },
                                        ad::wrt(p.b_Rsun, p.b_b0, p.b_z0, p.b_r0, p.b_psi0, p.b_psi1, p.b_xsi0), ad::at(x, y, z, p), out);
  return _filter_diff(_deriv);
}

#endif