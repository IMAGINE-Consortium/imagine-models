#ifndef SUN_H
#define SUN_H


#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

//Sun et al. A&A V.477 2008 ASS+RING model magnetic field


class SunMagneticField : public RegularVectorField  {
    protected:

    vector _at_position(const double &x, const double &y, const double &z, const SunMagneticField &p) const;

#if autodiff_FOUND
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, SunMagneticField &p) const;
#endif
    public:
        using RegularVectorField :: RegularVectorField;

        number b_Rsun = 8.5; 
        number b_R0 = 8.5;
        number b_B0 = 2.;
        number b_z0 = 1.;
        number b_Rc = 5.3;
        number b_Bc = 2.;
        number b_p = -12.;
        
        number bH_B0 = 2.;  // 10 in original publication, 2 in update https://arxiv.org/abs/1010.4394
        number bH_R0 = 4.;
        number bH_z0 = 1.5;
        number bH_z1a = 0.2;
        number bH_z1b = 0.4;

#if autodiff_FOUND
    const std::set<std::string> all_diff{"b_Rsun", "b_B0", "b_R0", "b_z0", "b_Rc", "b_Bc", "b_p", "bH_B0", "bH_R0", "bH_z0",  "bH_z1a", "bH_z1b"};
    std::set<std::string> active_diff{"b_Rsun", "b_B0", "b_R0", "b_z0", "b_Rc", "b_Bc", "b_p", "bH_B0", "bH_R0", "bH_z0",  "bH_z1a", "bH_z1b"};

    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
    {
        return _jac(x, y, z, *this);
    }
#endif
        vector at_position(const double &x, const double &y, const double &z) const {
            return _at_position(x, y, z, *this);
        }
 };

 #endif
