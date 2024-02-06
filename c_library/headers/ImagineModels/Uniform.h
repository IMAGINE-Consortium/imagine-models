#ifndef UNIFORM_H
#define UNIFORM_H

#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

class UniformMagneticField : public RegularVectorField
{
protected:
    vector _at_position(const double &x, const double &y, const double &z, const UniformMagneticField &p) const
    {
        return vector{{p.bx, p.by, p.bz}};
    }

#if autodiff_FOUND
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, UniformMagneticField &p) const
    {
        return Eigen::MatrixXd::Identity(3, 3);
    }
#endif
public:
    using RegularVectorField ::RegularVectorField;

    number bx = 0.;
    number by = 0.;
    number bz = 0.;
#if autodiff_FOUND
    const std::set<std::string> all_diff{"bx", "by", "bz"};
    std::set<std::string> active_diff{"bx", "by", "bz"};

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


class UniformDensityField : public RegularScalarField
{
protected:
    number _at_position(const double &x, const double &y, const double &z, const UniformDensityField &p) const
    {
        return p.n0;
    }

#if autodiff_FOUND
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, UniformDensityField &p) const
    {
        return Eigen::MatrixXd::Identity(1, 1);
    }
#endif
public:
    using RegularScalarField ::RegularScalarField;

    number n0 = 0.;

#if autodiff_FOUND
    const std::set<std::string> all_diff{"n0"};
    std::set<std::string> active_diff{"n0"};

    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
    {
        return _jac(x, y, z, *this);
    }
#endif
    number at_position(const double &x, const double &y, const double &z) const
    {
        return _at_position(x, y, z, *this);
    }
};

#endif