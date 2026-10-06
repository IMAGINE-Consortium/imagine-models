#ifndef TT_H
#define TT_H

#include <functional>
#include <cmath>

#include "ImagineModels/RegularField.h"

namespace imagine {

// Tinyakov and Tkachev (TT) https://arxiv.org/abs/astro-ph/0111305, implementation of https://arxiv.org/pdf/astro-ph/0510444.pdf (Kachelriess et al.)
class TTMagneticField : public RegularVectorField
{
protected:
    vector _at_position(const double &x, const double &y, const double &z, const TTMagneticField &p) const;

#if IMAGINE_HAS_AUTODIFF
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, TTMagneticField &p) const;
#endif
public:
    using RegularVectorField::RegularVectorField;

    double b_r_max = 20.; // kpc
    double b_r_min = 4.;  // kpc

    number b_Rsun = 8.5; // kpc
    number b_b0 = 1.4;   // muG
    number b_d = -0.5;   // kpc
    number b_z0 = 1.5;   // kpc, called h in original publication
    number b_p = -8;     // degree

#if IMAGINE_HAS_AUTODIFF
    const std::set<std::string> all_diff{"b_Rsun", "b_b0", "b_d", "b_z0", "b_p"};
    std::set<std::string> active_diff{"b_Rsun", "b_b0", "b_d", "b_z0", "b_p"};

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

}

#endif
