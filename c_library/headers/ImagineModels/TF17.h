#ifndef TERRALFERRIERE_H
#define TERRALFERRIERE_H

#include <functional>
#include <cmath>

#include "Field.h"
#include "RegularField.h"

// Terral, Ferriere 2017 - Constraints from Faraday rotation on the magnetic field structure in the galactic halo, DOI: 10.1051/0004-6361/201629572, arXiv:1611.10222, implementation adapted from CRPRopa

class TFMagneticField : public RegularVectorField
{
protected:
    vector _at_position(const double &x, const double &y, const double &z, const TFMagneticField &p) const;

#if autodiff_FOUND
    Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, TFMagneticField &p) const;
#endif
public:
    using RegularVectorField ::RegularVectorField;

    std::string activeDiskModel = "Ad1";
    const std::array<std::string, 3> possibleDiskModels{"Ad1", "Bd1", "Dd1"};
    std::string activeHaloModel = "C0";
    const std::array<std::string, 2> possibleHaloModels{"C0", "C1"};

    // disk parameters
    number a_disk = 0.9;         // kp**-2; not relevant for: Bd1, Dd1
    number z1_disk = 0;          // not relevant for: Ad1, Bd1
    number r1_disk = 3;          // kpc; // not relevant for: Dd1
    number B1_disk = 19.;        // muG;
    number L_disk = 0;           // not relevant for: Ad1, Bd1
    number phi_star_disk = -54.; // deg ;
    number H_disk = 0.0055;      // kpc; // not relevant for: Dd1

    // halo parameters
    number a_halo = 1.17;     // kp**-2;
    number z1_halo = 0.;      // kpc
    number B1_halo = 0.36;    // muG
    number L_halo = 3.0;      // kpc
    number phi_star_halo = 0; // deg

    // universal parameters
    number p_0 = -7.9; // deg;
    number H_p = 5.;   // kpc; Ad1
    number L_p = 50.;  // kpc

    // security to avoid 0 division
    double epsilon = 1e-16;

#if autodiff_FOUND
    const std::set<std::string> all_diff{"a_disk", "z1_disk", "r1_disk", "B1_disk", "L_disk", "phi_star_disk", "H_disk", "a_halo", "z1_halo", "B1_halo", "L_halo", "phi_star_halo", "p_0", "H_p", "L_p"};
    std::set<std::string> active_diff{"a_disk", "r1_disk", "B1_disk", "phi_star_disk", "H_disk", "a_halo", "z1_halo", "B1_halo", "L_halo", "p_0", "H_p", "L_p"};

    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
    {
        return _jac(x, y, z, *this);
    }
#endif
    vector at_position(const double &x, const double &y, const double &z) const
    {
        return _at_position(x, y, z, *this);
    }

    vector getDiskField(const double &r, const double &z, const double &phi, const double &sinPhi, const double &cosPhi, const TFMagneticField &p) const;

    vector getHaloField(const double &r, const double &z, const double &phi, const double &sinPhi, const double &cosPhi, const TFMagneticField &p) const;

    number azimuthalFieldComponent(const double &r, const double &z, const number &B_r, const number &B_z, const number &cp0, const TFMagneticField &p) const;

    number radialFieldScale(const number &B1, const number &phi_star, const number &z1, const double &phi, const double &r, const double &z, const number &cp0, const TFMagneticField &p) const;

    number shiftedWindingFunction(const number &r, const double &z, const number &cp0, const TFMagneticField &p) const;

    number zscale(const double &z, const TFMagneticField &p) const;

    void set_params(std::string dtype, std::string htype);
};

#endif