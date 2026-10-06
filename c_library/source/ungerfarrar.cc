/*
This file contains code adapted from Unger&Farrar 2024.

The original copyright statement is reproduced below:

BSD 2-Clause License

Copyright (c) 2024, Michael Unger and Glennys R. Farrar

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <limits>
#include <algorithm>
#include <stdexcept>

#include "ImagineModels/UngerFarrar.h"
#include "ImagineModels/helpers.h"
#include "ImagineModels/units.h"

namespace imagine {



void UFMagneticField::set_model(const std::string &model)
{
  if (std::find(available_models.begin(), available_models.end(), model) == available_models.end())
    throw std::invalid_argument("Unknown UF24 model '" + model + "'.");
  active_model = model;
  parameters = UFParameters<double>{};
  set_parameter_map(all_parameters.at(model));
}


template <typename T>
Vec3<T> UFMagneticField::field(const double &x, const double &y, const double &z, const UFParameters<T> &p) const
{
  Vec3<T> B_cart{{0., 0., 0.}};
  double squared_length = pow(x, 2) + pow(y, 2) + pow(z, 2);
  if (squared_length > pow(fMaxRadius, 2))
    return B_cart;
  else {
    const auto diskField = GetDiskField(x, y, z, p);
    const auto haloField = GetHaloField(x, y, z, p);
    for (size_t l = 0; l < 3; l++) {
      B_cart[l] = diskField[l] + haloField[l];
    }
    return B_cart;
  }
}

template <typename T>
Vec3<T> UFMagneticField::GetDiskField(const double &x, const double &y, const double &z, const UFParameters<T> &p) const
{
  if (active_model == "spur")
    return GetSpurField(x, y, z, p);
  else
    return GetSpiralField(x, y, z, p);
}


template <typename T>
Vec3<T> UFMagneticField::GetHaloField(const double &x, const double &y, const double &z, const UFParameters<T> &p) const
{
  if (active_model == "twistX")
    return GetTwistedHaloField(x, y, z, p);
  else {
    Vec3<T> B_cart_halo{{0., 0., 0.}};
    const auto poloidalHaloField = GetPoloidalHaloField(x, y, z, p);
    const auto toroidalHaloField = GetToroidalHaloField(x, y, z, p);
    for (size_t l = 0; l < 3; l++) {
      B_cart_halo[l] = toroidalHaloField[l] + poloidalHaloField[l];
    }
    return B_cart_halo;
  }
      
}


template <typename T>
Vec3<T> UFMagneticField::GetTwistedHaloField(const double x, const double y, const double z, const UFParameters<T> &p) const
{
  const double r = sqrt(x*x + y*y);
  const double cosPhi = r > std::numeric_limits<double>::min() ? x / r : 1;
  const double sinPhi = r > std::numeric_limits<double>::min() ? y / r : 0;

  Vec3<T> bXCart = GetPoloidalHaloField(x, y, z, p);
  Vec3<T> bXCartTmp{{bXCart[0], bXCart[1], bXCart[2]}};
  Vec3<T> bXCyl = Cart2Cyl(bXCartTmp, cosPhi, sinPhi);

  T bZ = bXCyl[2];
  T bR = bXCyl[0];

  T bPhi = 0;

  if (p.fTwistingTime != 0 && r != 0) {
    // radial rotation curve parameters (fit to Reid et al 2014)
    const double v0 = -240 * astro::kilometer/astro::second;
    const double r0 = 1.6; // kpc
    // vertical gradient (Levine+08)
    const double z0 = 10; //

    // Eq.(43)
    const double fr = 1 - exp(-r/r0);
    // Eq.(44)
    const double t0 = exp(2*abs(z)/z0);
    const double gz = 2 / (1 + t0);

    // Eq. (46)
    const double signZ = z < 0 ? -1 : 1;
    const double deltaZ =  -signZ * v0 * fr / z0  * t0 * pow(gz, 2);
    // Eq. (47)
    const double deltaR = v0 * ((1-fr)/r0 - fr/r) * gz;

    // Eq.(45)
    bPhi = (bZ * deltaZ + bR * deltaR) * p.fTwistingTime;

  }
  Vec3<T> bCylX{{bR, bPhi , bZ}};
  return Cyl2Cart<Vec3<T>>(bCylX, cosPhi, sinPhi);
}

template <typename T>
Vec3<T> UFMagneticField::GetToroidalHaloField(const double x, const double y, const double z, const UFParameters<T> &p) const
{
  const double r2 = x*x + y*y;
  const double r = sqrt(r2);
  const double absZ = abs(z);

  T b0 = z >= 0 ? p.fToroidalBN : p.fToroidalBS;
  T rh = p.fToroidalR;
  T z0 = p.fToroidalZ;
  T fwh = p.fToroidalW;
  //number sigmoidR = Sigmoid<number>(r, rh, fwh);
  T sigmoidR = 1 / (1 + exp(-(r-rh)/fwh));
  //number sigmoidZ = Sigmoid<number>(absZ, p.fDiskH, p.fDiskW);
  T sigmoidZ = 1 / (1 + exp(-(absZ-p.fDiskH)/p.fDiskW));

  // Eq. (21)
  T bPhi = b0 * (1. - sigmoidR) * sigmoidZ * exp(-absZ/z0);

  Vec3<T> bCyl{{0., bPhi, 0.}};
  const double cosPhi = r > std::numeric_limits<double>::min() ? x / r : 1;
  const double sinPhi = r > std::numeric_limits<double>::min() ? y / r : 0;
  return Cyl2Cart<Vec3<T>>(bCyl, cosPhi, sinPhi);
}

template <typename T>
Vec3<T> UFMagneticField::GetPoloidalHaloField(const double x, const double y, const double z, const UFParameters<T> &p) const
{
  const double r2 = x*x + y*y;
  const double r = std::sqrt(r2);

  T c = pow(p.fPoloidalA/p.fPoloidalZ, p.fPoloidalP);
  T a0p = pow(p.fPoloidalA, p.fPoloidalP);
  T rp = pow(r, p.fPoloidalP);
  T abszp = pow(abs(z), p.fPoloidalP);
  T cabszp = c*abszp;

  /*
    since $\sqrt{a^2 + b} - a$ is numerical unstable for $b\ll a$,
    we use $(\sqrt{a^2 + b} - a) \frac{\sqrt{a^2 + b} + a}{\sqrt{a^2
    + b} + a} = \frac{b}{\sqrt{a^2 + b} + a}$}
  */

  T t0 = a0p + cabszp - rp;
  T t1 = sqrt(pow(t0, 2) + 4*a0p*rp);
  T ap = 2*a0p*rp / (t1  + t0);

  T a = 0;
  if (ap < 0) {
    if (r > std::numeric_limits<double>::min()) {
      // this should never happen
      throw std::runtime_error("ap became negative and r is finite");
    }
    else
      a = 0;
  }
  else
    a = pow(ap, 1/p.fPoloidalP);

  // Eq.(29) and Eq.(32)
  T radialDependence =
    active_model == "expX" ?
    exp(-a/p.fPoloidalR) :
    //1 - Sigmoid<number>(a, p.fPoloidalR, p.fPoloidalW);
    1 - 1 / (1 + exp(-(a-p.fPoloidalR)/p.fPoloidalW));

  // Eq.(28)
  T Bzz = p.fPoloidalB * radialDependence;

  // (r/a)
  T rOverA =  1 / pow(2*a0p / (t1  + t0), 1/p.fPoloidalP);

  // Eq.(35) for p=n
  const double signZ = z < 0 ? -1 : 1;
  T Br =
    Bzz * c * a / rOverA * signZ * pow(abs(z), p.fPoloidalP - 1) / t1;

  // Eq.(36) for p=n
  T Bz = Bzz * pow(rOverA, p.fPoloidalP-2) * (ap + a0p) / t1;

  if (r < std::numeric_limits<double>::min())
    return Vec3<T>{{0., 0., Bz}};
  else {
    Vec3<T> bCylX{{Br, 0 , Bz}};
    const double cosPhi =  x / r;
    const double sinPhi =  y / r;
    return Cyl2Cart<Vec3<T>>(bCylX, cosPhi, sinPhi);
  }
}

template <typename T>
Vec3<T> UFMagneticField::GetSpurField(const double x, const double y, const double z, const UFParameters<T> &p) const
{
  // reference approximately at solar radius
  const double rRef = 8.2; //kpc

  auto fSinPitch = sin(p.fDiskPitch);
  auto fCosPitch = cos(p.fDiskPitch);
  auto fTanPitch = tan(p.fDiskPitch);
  // cylindrical coordinates
  const double r2 = x*x + y*y;
  const double r = sqrt(r2);
  if (r < std::numeric_limits<double>::min())
    return Vec3<T>{{0, 0, 0}};

  double phi = atan2(y, x);
  if (phi < 0)
    phi += num::twopi;

  T phiRef = p.fDiskPhase1;
  int iBest = -2;
  T bestDist = -1;
  for (int i = -1; i <= 1; ++i) {
    T pphi = phi - phiRef + i*num::twopi;
    T rr = rRef*exp(pphi * fTanPitch);
    if (bestDist < 0 || abs(r-rr) < bestDist) {
      bestDist =  abs(r-rr);
      iBest = i;
    }
  }
  if (iBest == 0) {
    T phi0 = phi - log(r/rRef) / fTanPitch;

    // Eq. (16)
    //number deltaPhi0 = DeltaPhi<number, number, number>(phiRef, phi0);
    T deltaPhi0 = acos(cos(phi0)*cos(phiRef) + sin(phi0)*sin(phiRef));
    T delta = deltaPhi0 / p.fSpurWidth;
    T B = p.fDiskB1 * exp(-0.5*pow(delta, 2));

    // Eq. (18)
    const double wS = 5*num::rad;
    T phiC = p.fSpurCenter;
    //number deltaPhiC = DeltaPhi<number, number, number>(phiC, phi);
    T deltaPhiC = acos(cos(phi)*cos(phiC) + sin(phi)*sin(phiC));
    T lC = p.fSpurLength;
    //number gS = 1 - Sigmoid<number>(abs(deltaPhiC), lC, wS);
    T gS = 1 - 1 / (1 + exp(-(abs(deltaPhiC)-lC)/wS));

    // Eq. (13)
    //number hd = 1 - Sigmoid<number>(abs(z), p.fDiskH, p.fDiskW);
    T hd = 1 - 1 / (1 + exp(-(abs(z)-p.fDiskH)/p.fDiskW));

    // Eq. (17)
    T bS = rRef/r * B * hd * gS;
    Vec3<T> bCyl{{bS * fSinPitch, bS * fCosPitch, 0.}};
    const double cosPhi = x / r;
    const double sinPhi = y / r;
    return Cyl2Cart<Vec3<T>>(bCyl, cosPhi, sinPhi);
  }
  else
    return Vec3<T>{{0, 0, 0}};

}

template <typename T>
Vec3<T> UFMagneticField::GetSpiralField(const double x, const double y, const double z, const UFParameters<T> &p) const
{
  // reference radius
  const double rRef = 5.; // kpc
  // inner boundary of spiral field
  const double rInner = 5; // kpc
  const double wInner = 0.5; // kpc
  // outer boundary of spiral field
  const double rOuter = 20; // kpc
  const double wOuter = 0.5; // kpc

  auto fSinPitch = sin(p.fDiskPitch);
  auto fCosPitch = cos(p.fDiskPitch);
  auto fTanPitch = tan(p.fDiskPitch);

  // cylindrical coordinates
  const double r2 = x*x + y*y;
  if (r2 == 0)
    return Vec3<T>{{0, 0, 0}};

  const double r = std::sqrt(r2);
  const double phi = std::atan2(y, x);

  // Eq.(13)
  //number hdz = 1 - Sigmoid(abs(z), fDiskH, fDiskW);
  T hdz = 1 - 1 / (1 + exp(-(abs(z)-p.fDiskH)/p.fDiskW));

  // Eq.(14) times rRef divided by r
  //const double rFacI = Sigmoid(r, rInner, wInner);
  const double rFacI = 1 / (1 + exp(-(r-rInner)/wInner));
  //const double rFacO = 1 - Sigmoid(r, rOuter, wOuter);
  const double rFacO = 1 - 1 / (1 + exp(-(r-rOuter)/wOuter));
  
  // (using lim r--> 0 (1-exp(-r^2))/r --> r - r^3/2 + ...)
  const double rFac =  r > 1e-5*astro::pc ? (1-exp(-r*r)) / r : r * (1 - r2/2);
  const double gdrTimesRrefByR = rRef * rFac * rFacO * rFacI;

  // Eq. (12)
  T phi0 = phi - log(r/rRef) / fTanPitch;

  // Eq. (10)
  T b =
    p.fDiskB1 * cos(1 * (phi0 - p.fDiskPhase1)) +
    p.fDiskB2 * cos(2 * (phi0 - p.fDiskPhase2)) +
    p.fDiskB3 * cos(3 * (phi0 - p.fDiskPhase3));

  // Eq. (11)
  T fac = hdz * gdrTimesRrefByR;
  Vec3<T> bCyl{{ b * fac * fSinPitch,
      b * fac * fCosPitch,
      0.}};

  const double cosPhi = x / r;
  const double sinPhi = y / r;
  return Cyl2Cart<Vec3<T>>(bCyl, cosPhi, sinPhi);
}


IMAGINE_INSTANTIATE_VECTOR_MODEL(UFMagneticField)

}
