// Reference: Kachelriess et al. 2007, arXiv:astro-ph/0510444 (Sec. II C), after Prouza & Smida 2003, arXiv:astro-ph/0307165
// Deviations:
// - version of Kachelriess et al. 2007, not the original of Prouza & Smida 2003
// - disk amplitude constant for r < 4 kpc as in the TT model (the paper leaves the inner disk open)
// - r_max = 20 kpc cut applied to the disk only
// - dipole core (R < 0.5 kpc): B = (0, 0, -100 muG)
// - the halo's solar-circle radius h_R0 is a setting

#pragma once

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define PS_PARAMETERS(X)     \
    X(b_Rsun, 8.5)           \
    X(b_b0, 2.)              \
    X(b_d, -0.5)             \
    X(b_z0, 0.2)             \
    X(b_p, -8.) /* degree */ \
    X(h_b0, 1.5)             \
    X(h_z0, 1.5)             \
    X(h_w, 0.3)              \
    X(d_mu, 123.) /* muG kpc^3 */

IMAGINE_PARAMETERS(PSParameters, PS_PARAMETERS)

class PSMagneticField : public RegularVectorModel<PSMagneticField, PSParameters> {
public:
    bool do_disk = true;
    bool do_halo = true;
    bool do_dipole = true;

    double b_r_max = 20.;   // kpc
    double b_r_min = 4.;    // kpc
    double h_R0 = 8.5;      // kpc
    double d_r_core = 0.5;  // kpc
    double d_b_core = -100; // muG

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const PSParameters<T> &p) const;
};

}
