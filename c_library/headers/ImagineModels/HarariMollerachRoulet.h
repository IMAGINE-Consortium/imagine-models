#ifndef HMR_H
#define HMR_H

#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

// Harari, Mollerach, Roulet (HMR) see https://arxiv.org/abs/astro-ph/9906309, implementation of https://arxiv.org/pdf/astro-ph/0510444.pdf

class HMRMagneticField : public RegularVectorField
{
protected:
    vector _at_position(const double &x, const double &y, const double &z, const HMRMagneticField &p) const;

#if autodiff_FOUND
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, HMRMagneticField &p) const;
#endif
public:
    using RegularVectorField ::RegularVectorField;

    double b_r_max = 20.; // kpc

    number b_Rsun = 8.5;       // kpc
    number b_z1 = 0.3;         // kpc
    number b_z2 = 4.;          // kpc
    number b_r1 = 2.;          // kpc
    number b_p = -10;          // degree
    number b_epsilon0 = 10.55; // kpc

#if autodiff_FOUND
    const std::set<std::string> all_diff{"b_Rsun", "b_z1", "b_z2", "b_r1", "b_p", "b_epsilon0"};
    std::set<std::string> active_diff{"b_Rsun", "b_z1", "b_z2", "b_r1", "b_p", "b_epsilon0"};
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