#ifndef UNIFORM_H
#define UNIFORM_H

#include "ImagineModels/RegularModel.h"

namespace imagine {

#define UNIFORM_MAGNETIC_PARAMETERS(X) \
    X(bx, 0.)                          \
    X(by, 0.)                          \
    X(bz, 0.)

IMAGINE_PARAMETERS(UniformMagneticParameters, UNIFORM_MAGNETIC_PARAMETERS)

class UniformMagneticField : public RegularVectorModel<UniformMagneticField, UniformMagneticParameters>
{
public:
    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const UniformMagneticParameters<T> &p) const
    {
        return {p.bx, p.by, p.bz};
    }
};

#define UNIFORM_DENSITY_PARAMETERS(X) \
    X(n0, 0.)

IMAGINE_PARAMETERS(UniformDensityParameters, UNIFORM_DENSITY_PARAMETERS)

class UniformDensityField : public RegularScalarModel<UniformDensityField, UniformDensityParameters>
{
public:
    template <typename T>
    T field(const double &x, const double &y, const double &z, const UniformDensityParameters<T> &p) const
    {
        return p.n0;
    }
};

}

#endif
