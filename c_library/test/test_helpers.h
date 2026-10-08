#ifndef IMAGINE_TEST_HELPERS_H
#define IMAGINE_TEST_HELPERS_H

#include <array>
#include <cmath>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "ImagineModels/RegularModels.h"

namespace imagine::test {

using VectorModels =
    std::tuple<ArchimedeanMagneticField, FauvetMagneticField, HanMagneticField, HelixMagneticField, HMRMagneticField,
               JaffeMagneticField, JF12MagneticField, PshirkovMagneticField, StanevBSSMagneticField, SunMagneticField,
               SVT22MagneticField, TFMagneticField, TTMagneticField, UFMagneticField, UniformMagneticField,
               WMAPMagneticField, XH24MagneticField>;
using ScalarModels = std::tuple<UniformDensityField, YMW16>;
using AllModels = decltype(std::tuple_cat(std::declval<VectorModels>(), std::declval<ScalarModels>()));

using Position = std::array<double, 3>;

inline const std::vector<Position> positions = {
    {0., 0., 0.},        {-8.5, 0., 0.},     {8.5, 0., 0.},    {0., 8.5, 0.},  {-8.5, 0., 1.},    {-8.5, 0., -1.},
    {3., 4., 0.2},       {-1., -1., -1.},    {12., -9., 2.},   {0., 0., 5.},   {-4.2, 6.1, -0.3}, {15.3, 2.2, 3.9},
    {-17.8, -12.4, 1.1}, {5.5, -14.9, -4.2}, {1.2, 0.7, 0.05}, {-25., 3., 0.5}};

inline const std::vector<Position> z_axis = {{0., 0., -3.}, {0., 0., -.5}, {0., 0., 0.}, {0., 0., .5}, {0., 0., 3.}};

inline std::vector<double> as_values(const Vec3<double> &v) {
    return {v.begin(), v.end()};
}
inline std::vector<double> as_values(double v) {
    return {v};
}

template <typename Model> std::vector<double> value_at(const Model &model, const Position &p) {
    return as_values(model.at_position(p[0], p[1], p[2]));
}

inline bool all_finite(const std::vector<double> &values) {
    for (double v : values)
        if (!std::isfinite(v))
            return false;
    return true;
}

inline std::string to_string(const Position &p) {
    return "(" + std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " + std::to_string(p[2]) + ")";
}

}

#endif
