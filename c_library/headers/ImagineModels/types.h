#ifndef IMAGINE_TYPES_H
#define IMAGINE_TYPES_H

#include <array>
#include <cmath>
#include <cstdlib>

#include "ImagineModels/config.h"

#if IMAGINE_HAS_AUTODIFF
    #include <autodiff/forward/real.hpp>
    #include <autodiff/forward/real/eigen.hpp>
#endif

namespace imagine {

using std::abs;

#if IMAGINE_HAS_AUTODIFF
    namespace ad = autodiff;
#endif

template <typename T>
using Vec3 = std::array<T, 3>;

template <typename T>
using Scalar = T;

}

#endif
