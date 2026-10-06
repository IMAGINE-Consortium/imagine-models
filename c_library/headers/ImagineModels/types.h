#ifndef IMAGINE_TYPES_H
#define IMAGINE_TYPES_H

#include <array>

#include "ImagineModels/config.h"

#if IMAGINE_HAS_AUTODIFF
    #include <autodiff/forward/real.hpp>
    #include <autodiff/forward/dual.hpp>
    #include <autodiff/forward/real/eigen.hpp>
#endif

namespace imagine {

#if IMAGINE_HAS_AUTODIFF
    namespace ad = autodiff;
    typedef ad::real number;
    typedef ad::VectorXreal vector;
#else
    typedef double number;  // only used for differentiable numbers! 
    typedef std::array<double, 3> vector;
#endif

}

#endif
