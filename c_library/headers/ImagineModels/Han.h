#ifndef HAN_H
#define HAN_H


#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

//J. L. Han et al 2018 ApJS 234 11


class HanMagneticField : public RegularVectorField  {
    protected:

    vector _at_position(const double &x, const double &y, const double &z, const HanMagneticField &p) const;

#if autodiff_FOUND
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, HanMagneticField &p) const;
#endif
    public:
        using RegularVectorField :: RegularVectorField;

        number B_p = 11; // pitch angle 
        number A = 5.;
        number H = 0.4;
        number B_s1 = 4.5; 
        number B_s2 = -3.0;
        number B_s3 = 6.3;
        number B_s4 = -4.7;
        number B_s5 = 3.3;
        number B_s6 = -8.7;

        double R_min = 3.;
        double R_max = 15.;
        std::array<double, 7> R_s{3.0, 4.1, 4.9, 6.1, 7.5, 8.5, 10.5};
#if autodiff_FOUND
    const std::set<std::string> all_diff{"B_p", "A", "H", "B_s1", "B_s2", "B_s3", "B_s4", "B_s5", "B_s6"};
    std::set<std::string> active_diff{"B_p", "A", "H", "B_s1", "B_s2", "B_s3", "B_s4", "B_s5", "B_s6"};

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
