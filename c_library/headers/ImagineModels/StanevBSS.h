#ifndef STANEVBSS_H
#define STANEVBSS_H


#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"


// StanevBSS (HMR) see https://arxiv.org/abs/astro-ph/9607086

class StanevBSSMagneticField : public RegularVectorField  {
    protected:
protected:
    vector _at_position(const double &x, const double &y, const double &z, const StanevBSSMagneticField &p) const;

#if autodiff_FOUND
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, StanevBSSMagneticField &p) const;
#endif
    public:
        using RegularVectorField :: RegularVectorField;

        double b_r_max = 20.; // kpc
        double b_r_min = 4.; // kpc

        number b_z01 = 1.; // kpc
        number b_z02 = 4.; // kpc
        number b_z0_border = 0.5; // kpc
        number b_r0 = 10.55; // kpc
        number b_p = -10; // degree
        number b_Rsun = 8.5; // kpc

        number b_phi0 = M_PI; // radians
    
#if autodiff_FOUND
    const std::set<std::string> all_diff{"b_Rsun", "b_b0", "b_z0", "b_r0", "b_p", "b_phi0"};
    std::set<std::string> active_diff{"b_Rsun", "b_b0", "b_z0", "b_r0", "b_p", "b_phi0"};

    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
    {
        return _jac(x, y, z, *this);
    }
#endif

    vector at_position(const double &x, const double &y, const double &z) const
    {
        return _at_position(x, y, z, *this);
    }
 };

 #endif