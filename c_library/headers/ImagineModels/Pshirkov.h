#ifndef PSHIRKOV_H
#define PSHIRKOV_H

#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

class PshirkovMagneticField : public RegularVectorField
{
protected:
    vector _at_position(const double &xx, const double &yy, const double &zz, const PshirkovMagneticField &p) const;

#if autodiff_FOUND
    Eigen::MatrixXd _jac(const double &xx, const double &yy, const double &zz, PshirkovMagneticField &p) const;
#endif

public:
    using RegularVectorField ::RegularVectorField;

	bool useASS = false;  // switch for axisymmetric spiral field (ASS)
	bool useBSS = true;  // switch for bisymmetric spiral field (BSS)
	bool useHalo = true; // switch for halo field

    // differentiable parameters

	// disk parameters
	number pitch = -6;  //pitch angle parameters, deg (paper usues -5 for ASS and -6 for BSS)
	number d = - 0.6;     // distance to first field reversal, kpc
	number R_sun = 8.5; // distance between sun and galactic center, kpc
	double R_c = 5.0;   // radius of central region, kpc
	number z0_D = 1.0;    // vertical thickness in the galactic disk, kpc
	number B0_D = 2.0;    // magnetic field scale, muG

	// halo parameters
	number z0_H = 1.3;  // halo vertical position, kpc
	number R0_H = 8.0;  // halo radial position, kpc
	number B0_Hn = 4.0; // halo magnetic field scale (north), muG
	number B0_Hs = 4.0; // halo magnetic field scale (south), muG (paper usues 2 for ASS and 4 for BSS)
	number z11_H = 0.25; // halo vertical thickness towards disc, kpc
	number z12_H = 0.4; // halo vertical thickness off the disk, kpc


    #if autodiff_FOUND
        const std::set<std::string> all_diff{"pitch", "d", "R_sun", "z0_D", "B0_D", "z0_H", "R0_H", "B0_Hn", "B0_Hs", "z11_H", "z12_H"};
        std::set<std::string> active_diff{"pitch", "d", "R_sun","z0_D", "B0_D", "z0_H", "R0_H", "B0_Hn", "B0_Hs", "z11_H", "z12_H"};
    #endif

    vector at_position(const double &x, const double &y, const double &z) const
    {
        return _at_position(x, y, z, *this);
    }

#if autodiff_FOUND
    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
    {
        return _jac(x, y, z, *this);
    }
#endif
};

#endif
