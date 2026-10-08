#pragma once

#include <array>
#include <cmath>
#include <cstdlib>

#include "ImagineModels/config.h"

#if IMAGINE_HAS_AUTODIFF
#include <Eigen/Core>
#include <autodiff/forward/real.hpp>
#endif

namespace imagine {

using std::abs;

#if IMAGINE_HAS_AUTODIFF
namespace ad = autodiff;
#endif

template <typename T> using Vec3 = std::array<T, 3>;

template <typename T> using Scalar = T;

inline double value(double x) {
    return x;
}

#if IMAGINE_HAS_AUTODIFF
inline double value(const ad::real &x) {
    return x[0];
}
#endif

}
