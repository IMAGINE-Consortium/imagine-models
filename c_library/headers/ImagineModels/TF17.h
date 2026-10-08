// Reference: Terral & Ferriere 2017, arXiv:1611.10222
// Based on: CRPropa (TF17Field)
// Deviations:
// - lower limits of Table 2 used as values, as in CRPropa

#pragma once

#include <cmath>
#include <functional>

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define TF17_PARAMETERS(X)                        \
    X(a_disk, 0.9)         /* kpc^-2; Ad1 only */ \
    X(z1_disk, 0)          /* Dd1 only */         \
    X(r1_disk, 3)          /* kpc; Ad1, Bd1 */    \
    X(B1_disk, 19.)        /* muG */              \
    X(L_disk, 0)           /* Dd1 only */         \
    X(phi_star_disk, -54.) /* deg */              \
    X(H_disk, 0.055)       /* kpc; Ad1, Bd1 */    \
    X(a_halo, 1.17)        /* kpc^-2 */           \
    X(z1_halo, 0.)         /* kpc */              \
    X(B1_halo, 0.36)       /* muG */              \
    X(L_halo, 3.0)         /* kpc */              \
    X(phi_star_halo, 0)    /* deg */              \
    X(p_0, -7.9)           /* deg */              \
    X(H_p, 5.)             /* kpc; Ad1 */         \
    X(L_p, 50.)            /* kpc */

IMAGINE_PARAMETERS(TF17Parameters, TF17_PARAMETERS)

class TF17MagneticField : public RegularVectorModel<TF17MagneticField, TF17Parameters> {
public:
    const std::array<std::string, 3> available_disk_models{"Ad1", "Bd1", "Dd1"};
    const std::array<std::string, 2> available_halo_models{"C0", "C1"};

    explicit TF17MagneticField(const std::string &disk_model = "Ad1", const std::string &halo_model = "C0") {
        set_model(disk_model, halo_model);
    }

    void set_model(const std::string &disk_model, const std::string &halo_model);
    const std::string &disk_model() const { return active_disk_model; }
    const std::string &halo_model() const { return active_halo_model; }

    // avoids division by zero
    double epsilon = 1e-16;

    template <typename T>
    Vec3<T> getDiskField(const double &r, const double &z, const double &phi, const double &sinPhi,
                         const double &cosPhi, const TF17Parameters<T> &p) const;

    template <typename T>
    Vec3<T> getHaloField(const double &r, const double &z, const double &phi, const double &sinPhi,
                         const double &cosPhi, const TF17Parameters<T> &p) const;

    template <typename T>
    T azimuthalFieldComponent(const double &r, const double &z, const T &B_r, const T &B_z, const T &cp0,
                              const TF17Parameters<T> &p) const;

    template <typename T>
    T radialFieldScale(const T &B1, const T &phi_star, const T &z1, const double &phi, const double &r, const double &z,
                       const T &cp0, const TF17Parameters<T> &p) const;

    template <typename T>
    T shiftedWindingFunction(const T &r, const double &z, const T &cp0, const TF17Parameters<T> &p) const;

    template <typename T> T zscale(const double &z, const TF17Parameters<T> &p) const;

    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const TF17Parameters<T> &p) const;

private:
    std::string active_disk_model = "Ad1";
    std::string active_halo_model = "C0";
};

}
