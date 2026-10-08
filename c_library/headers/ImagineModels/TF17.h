#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

// Terral, Ferriere 2017 - Constraints from Faraday rotation on the magnetic field structure in the galactic halo,
// DOI: 10.1051/0004-6361/201629572, arXiv:1611.10222, implementation adapted from CRPRopa

#define TF17_PARAMETERS(X)                                          \
    X(a_disk, 0.9)         /* kp**-2; not relevant for: Bd1, Dd1 */ \
    X(z1_disk, 0)          /* not relevant for: Ad1, Bd1 */         \
    X(r1_disk, 3)          /* kpc; // not relevant for: Dd1 */      \
    X(B1_disk, 19.)        /* muG; */                               \
    X(L_disk, 0)           /* not relevant for: Ad1, Bd1 */         \
    X(phi_star_disk, -54.) /* deg ; */                              \
    X(H_disk, 0.055)       /* kpc; // not relevant for: Dd1 */      \
    X(a_halo, 1.17)        /* kp**-2; */                            \
    X(z1_halo, 0.)         /* kpc */                                \
    X(B1_halo, 0.36)       /* muG */                                \
    X(L_halo, 3.0)         /* kpc */                                \
    X(phi_star_halo, 0)    /* deg */                                \
    X(p_0, -7.9)           /* deg; */                               \
    X(H_p, 5.)             /* kpc; Ad1 */                           \
    X(L_p, 50.)            /* kpc */

IMAGINE_PARAMETERS(TFParameters, TF17_PARAMETERS)

class TFMagneticField : public RegularVectorModel<TFMagneticField, TFParameters> {
public:
    const std::array<std::string, 3> available_disk_models{"Ad1", "Bd1", "Dd1"};
    const std::array<std::string, 2> available_halo_models{"C0", "C1"};

    explicit TFMagneticField(const std::string &disk_model = "Ad1", const std::string &halo_model = "C0") {
        set_model(disk_model, halo_model);
    }

    void set_model(const std::string &disk_model, const std::string &halo_model);
    const std::string &disk_model() const { return active_disk_model; }
    const std::string &halo_model() const { return active_halo_model; }

    // security to avoid 0 division
    double epsilon = 1e-16;

    template <typename T>
    Vec3<T> getDiskField(const double &r, const double &z, const double &phi, const double &sinPhi,
                         const double &cosPhi, const TFParameters<T> &p) const;

    template <typename T>
    Vec3<T> getHaloField(const double &r, const double &z, const double &phi, const double &sinPhi,
                         const double &cosPhi, const TFParameters<T> &p) const;

    template <typename T>
    T azimuthalFieldComponent(const double &r, const double &z, const T &B_r, const T &B_z, const T &cp0,
                              const TFParameters<T> &p) const;

    template <typename T>
    T radialFieldScale(const T &B1, const T &phi_star, const T &z1, const double &phi, const double &r, const double &z,
                       const T &cp0, const TFParameters<T> &p) const;

    template <typename T>
    T shiftedWindingFunction(const T &r, const double &z, const T &cp0, const TFParameters<T> &p) const;

    template <typename T> T zscale(const double &z, const TFParameters<T> &p) const;

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const TFParameters<T> &p) const;

private:
    std::string active_disk_model = "Ad1";
    std::string active_halo_model = "C0";
};

}
