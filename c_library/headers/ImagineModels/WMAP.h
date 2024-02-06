#ifndef WMAP_H
#define WMAP_H

#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

// WMAP magnetic field

class WMAPMagneticField : public RegularVectorField
{
protected:
    vector _at_position(const double &x, const double &y, const double &z, const WMAPMagneticField &p) const;

#if autodiff_FOUND
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, WMAPMagneticField &p) const;
#endif
public:
    using RegularVectorField ::RegularVectorField;

    double b_r_max = 20.; // kpc
    double b_r_min = 3.;  // kpc

    number b_Rsun = 8.; // kpc
    number b_b0 = 6.;    // muG  -> not given in original paper? Could also be 3 according to https://www.aanda.org/articles/aa/full_html/2010/14/aa12733-09/aa12733-09.html
    number b_z0 = 1.;    // kpc
    number b_r0 = 8.;    // kpc
    number b_psi0 = 35;  // degree
    number b_psi1 = 0.9; // degree
    number b_xsi0 = 25;  // degree

    bool anti = false;

#if autodiff_FOUND
    const std::set<std::string> all_diff{"b_Rsun", "b_b0", "b_z0", "b_r0", "b_psi0", "b_psi1", "b_xsi0"};
    std::set<std::string> active_diff{"b_Rsun", "b_b0", "b_z0", "b_r0", "b_psi0", "b_psi1", "b_xsi0"};

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