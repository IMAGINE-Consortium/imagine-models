#pragma once

namespace imagine {
namespace units {

constexpr double pi = 3.141592653589793238462643383279502884197;
constexpr double twopi = 2. * pi;
constexpr double halfpi = pi / 2.;
constexpr double deg = pi / 180.;

constexpr double kpc = 1.;
constexpr double pc = 1e-3 * kpc;
constexpr double Gpc = 1e6 * kpc;
constexpr double microgauss = 1.;
constexpr double megayear = 1.;
constexpr double second = megayear / (1e6 * 60 * 60 * 24 * 365.25);
constexpr double kilometer = kpc / 3.0856775807e+16;

}
}
