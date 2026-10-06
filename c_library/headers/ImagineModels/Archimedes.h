#ifndef ARCHIMEDES_H
#define ARCHIMEDES_H


#include <functional>
#include <cmath>

#include "ImagineModels/Field.h"
#include "ImagineModels/RegularField.h"

namespace imagine {

//simple archimdeean sprial, implementation based on CRPropa


class ArchimedeanMagneticField : public RegularVectorField  {
    protected:

    vector _at_position(const double &x, const double &y, const double &z, const ArchimedeanMagneticField &p) const;

#if IMAGINE_HAS_AUTODIFF
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, ArchimedeanMagneticField &p) const;
#endif
    public:
        using RegularVectorField :: RegularVectorField;

        number R_0 = 3; 
        number Omega = 1.;
        number v_w = 0.4;
        number B_0 = 1.;


#if IMAGINE_HAS_AUTODIFF
    const std::set<std::string> all_diff{"R_0", "Omega", "v_w", "B_0"};
    std::set<std::string> active_diff{"R_0", "Omega", "v_w", "B_0"};

    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
    {
        return _jac(x, y, z, *this);
    }
#endif
        vector at_position(const double &x, const double &y, const double &z) const {
            return _at_position(x, y, z, *this);
        }
 };

}

#endif
