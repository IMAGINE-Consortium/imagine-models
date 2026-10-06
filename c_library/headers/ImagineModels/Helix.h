#ifndef HELIX_H
#define HELIX_H

#include <functional>
#include <cmath>

#include "ImagineModels/RegularField.h"

namespace imagine {

class HelixMagneticField : public RegularVectorField
{
protected:
    vector _at_position(const double &xx, const double &yy, const double &zz, const HelixMagneticField &p) const;

#if IMAGINE_HAS_AUTODIFF
    Eigen::MatrixXd _jac(const double &xx, const double &yy, const double &zz, HelixMagneticField &p) const;
#endif

public:
    using RegularVectorField ::RegularVectorField;

    // differentiable parameters
    number ampx = 1.;
    number ampy = 1.;
    number ampz = 1.;

    #if IMAGINE_HAS_AUTODIFF
        const std::set<std::string> all_diff{"ampx", "ampy", "ampz"};
        std::set<std::string> active_diff{"ampx", "ampy", "ampz"};
    #endif

    // non_differentiable parameters
    double rmax = 20.;
    double rmin = 1.;

    vector at_position(const double &x, const double &y, const double &z) const
    {
        return _at_position(x, y, z, *this);
    }

#if IMAGINE_HAS_AUTODIFF
    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
    {
        return _jac(x, y, z, *this);
    }
#endif
};

/*THIS IS A ALTERNATIVE IMPLEMENTATION, WHICH COMPILES, BUT WHICH I DONT KNOW HOW TO BIND TO PYTHON IF WE WANT TO BE ABLE TO CHANGE THE PARAMETERS WRT DIFFERENTIATE AT RUNTIME. THE PROBLEM IS THAT AUTODIFF WANTS TO KNOW THIS AT COMPILE TIME. IF IT WOULD WORK, THE DERIVATIVE OWULD LIKELY BE FASTER, SINCE WE WOULDN'T NEED TO "FILTER".
            template<typename... Params>
            Eigen::MatrixXd _derivative_alt(const double &x, const double &y, const double &z,  HelixMagneticFields &p, Params... parameters) {
                vector out;

                Eigen::MatrixXd deriv = ad::jacobian([&](auto _x, auto _y, auto _z, auto _p) {return this->_at_position(_x, _y, _z, _p);},
                    ad::wrt(parameters...), ad::at(x, y, z, p), out);
                    return deriv;
                }

                        template<typename... Params>
            Eigen::MatrixXd derivative(const double &x, const double &y, const double &z, Params... parameters) {
                return _derivative_alt(x, y, z, this->param, parameters...);
            }

*/

}

#endif
