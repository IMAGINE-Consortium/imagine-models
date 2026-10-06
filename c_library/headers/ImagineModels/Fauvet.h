#ifndef FAUVET_H
#define FAUVET_H

#include <functional>
#include <cmath>

#include "ImagineModels/Field.h"
#include "ImagineModels/RegularField.h"

namespace imagine {

// Fauvet magnetic field

class FauvetMagneticField : public RegularVectorField
{
protected:
    vector _at_position(const double &x, const double &y, const double &z, const FauvetMagneticField &p) const;

#if IMAGINE_HAS_AUTODIFF
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, FauvetMagneticField &p) const;
#endif

public:
    using RegularVectorField ::RegularVectorField;

    double b_r_max = 20.; // kpc
    double b_r_min = 3.;  // kpc

    number b_b0 = 7.1;    // muG
    number b_z0 = 1.;     // kpc
    number b_r0 = 8.;     // kpc
    number b_p = -26.1;   // degree
    number b_chi0 = 22.4; // degree

    number h_b0 = 1.;  // muG
    number h_z0 = 1.5; // kpc
    number h_r0 = 4.;  // kpc
    number h_z1a = .2; // kpc
    number h_z1b = .4; // kpc

#if IMAGINE_HAS_AUTODIFF
    const std::set<std::string> all_diff{"b_b0", "b_z0", "b_r0", "b_p", "b_chi0", "h_b0", "h_z0", "h_r0", "h_z1a", "h_z1b"};
    std::set<std::string> active_diff{"b_b0", "b_z0", "b_r0", "b_p", "b_chi0", "h_b0", "h_z0", "h_r0", "h_z1a", "h_z1b"};

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