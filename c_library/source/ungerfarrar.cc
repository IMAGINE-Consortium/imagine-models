/*
This file contains code adapted from 

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

#include "UngerFarrar.h"
#include "helpers.h"

vector UFMagneticField::_at_position(const double &x, const double &y, const double &z, const UFMagneticField &p) const
{
  vector B_cart{{0., 0., 0.}};
  double squared_length = x**2 + y**2 + z**2;
  if (squared_length > p.fMaxRadiusSquared)
    return Vector3(0, 0, 0);
  else {
    const auto diskField = GetDiskField(x, y, z);
    const auto haloField = GetHaloField(x, y, z);
    return (diskField + haloField);
  }
}

vector UFMagneticField::GetDiskField(const double &x, const double &y, const double &z, const UFMagneticField &p)
  const
{
  if (p.fModelType == spur)
    return GetSpurField(x, y, z, const UFMagneticField &p);
  else
    return GetSpiralField(x, y, z, const UFMagneticField &p);
}


vector UFMagneticField::GetHaloField(const double &x, const double &y, const double &z, const UFMagneticField &p)
  const
{
  if (fModelType == twistX)
    return GetTwistedHaloField(x, y, z);
  else
    return
      GetToroidalHaloField(x, y, z) +
      GetPoloidalHaloField(x, y, z);
}


vector UFMagneticField::GetTwistedHaloField(const double x, const double y, const double z, const UFMagneticField &p)
  const
{
  const double r = sqrt(x*x + y*y);
  const double cosPhi = r > std::numeric_limits<double>::min() ? x / r : 1;
  const double sinPhi = r > std::numeric_limits<double>::min() ? y / r : 0;

  const Vector3 bXCart = GetPoloidalHaloField(x, y, z);
  const double bXCartTmp[3] = {bXCart.x, bXCart.y, bXCart.z};
  const Vector3 bXCyl = utl::CartToCyl(bXCartTmp, cosPhi, sinPhi);

  const double bZ = bXCyl.z;
  const double bR = bXCyl.x;

  double bPhi = 0;

  if (fTwistingTime != 0 && r != 0) {
    // radial rotation curve parameters (fit to Reid et al 2014)
    const double v0 = -240 * utl::kilometer/utl::second;
    const double r0 = 1.6 * utl::kpc;
    // vertical gradient (Levine+08)
    const double z0 = 10 * utl::kpc;

    // Eq.(43)
    const double fr = 1 - exp(-r/r0);
    // Eq.(44)
    const double t0 = exp(2*std::abs(z)/z0);
    const double gz = 2 / (1 + t0);

    // Eq. (46)
    const double signZ = z < 0 ? -1 : 1;
    const double deltaZ =  -signZ * v0 * fr / z0  * t0 * pow(gz, 2);
    // Eq. (47)
    const double deltaR = v0 * ((1-fr)/r0 - fr/r) * gz;

    // Eq.(45)
    bPhi = (bZ * deltaZ + bR * deltaR) * fTwistingTime;

  }
  const double bCylX[3] = {bR, bPhi , bZ};
  return utl::CylToCart(bCylX, cosPhi, sinPhi);
}

vector UFMagneticField::GetToroidalHaloField(const double x, const double y, const double z, const UFMagneticField &p)
  const
{
  const double r2 = x*x + y*y;
  const double r = sqrt(r2);
  const double absZ = std::abs(z);

  const double b0 = z >= 0 ? fToroidalBN : fToroidalBS;
  const double rh = fToroidalR;
  const double z0 = fToroidalZ;
  const double fwh = fToroidalW;
  const double sigmoidR = utl::Sigmoid(r, rh, fwh);
  const double sigmoidZ = utl::Sigmoid(absZ, fDiskH, fDiskW);

  // Eq. (21)
  const double bPhi = b0 * (1. - sigmoidR) * sigmoidZ * exp(-absZ/z0);

  const double bCyl[3] = {0, bPhi, 0};
  const double cosPhi = r > std::numeric_limits<double>::min() ? x / r : 1;
  const double sinPhi = r > std::numeric_limits<double>::min() ? y / r : 0;
  return utl::CylToCart(bCyl, cosPhi, sinPhi);
}

vector UFMagneticField::GetPoloidalHaloField(const double x, const double y, const double z, const UFMagneticField &p)
  const
{
  const double r2 = x*x + y*y;
  const double r = sqrt(r2);

  const double c = pow(fPoloidalA/fPoloidalZ, fPoloidalP);
  const double a0p = pow(fPoloidalA, fPoloidalP);
  const double rp = pow(r, fPoloidalP);
  const double abszp = pow(std::abs(z), fPoloidalP);
  const double cabszp = c*abszp;

  /*
    since $\sqrt{a^2 + b} - a$ is numerical unstable for $b\ll a$,
    we use $(\sqrt{a^2 + b} - a) \frac{\sqrt{a^2 + b} + a}{\sqrt{a^2
    + b} + a} = \frac{b}{\sqrt{a^2 + b} + a}$}
  */

  const double t0 = a0p + cabszp - rp;
  const double t1 = sqrt(pow(t0, 2) + 4*a0p*rp);
  const double ap = 2*a0p*rp / (t1  + t0);

  double a = 0;
  if (ap < 0) {
    if (r > std::numeric_limits<double>::min()) {
      // this should never happen
      throw std::runtime_error("ap = " + std::to_string(ap));
    }
    else
      a = 0;
  }
  else
    a = pow(ap, 1/fPoloidalP);

  // Eq.(29) and Eq.(32)
  const double radialDependence =
    fModelType == expX ?
    exp(-a/fPoloidalR) :
    1 - utl::Sigmoid(a, fPoloidalR, fPoloidalW);

  // Eq.(28)
  const double Bzz = fPoloidalB * radialDependence;

  // (r/a)
  const double rOverA =  1 / pow(2*a0p / (t1  + t0), 1/fPoloidalP);

  // Eq.(35) for p=n
  const double signZ = z < 0 ? -1 : 1;
  const double Br =
    Bzz * c * a / rOverA * signZ * pow(std::abs(z), fPoloidalP - 1) / t1;

  // Eq.(36) for p=n
  const double Bz = Bzz * pow(rOverA, fPoloidalP-2) * (ap + a0p) / t1;

  if (r < std::numeric_limits<double>::min())
    return Vector3(0, 0, Bz);
  else {
    const double bCylX[3] = {Br, 0 , Bz};
    const double cosPhi =  x / r;
    const double sinPhi =  y / r;
    return utl::CylToCart(bCylX, cosPhi, sinPhi);
  }
}

vector UFMagneticField::GetSpurField(const double x, const double y, const double z, const UFMagneticField &p)
  const
{
  // reference approximately at solar radius
  const double rRef = 8.2*utl::kpc;

  // cylindrical coordinates
  const double r2 = x*x + y*y;
  const double r = sqrt(r2);
  if (r < std::numeric_limits<double>::min())
    return Vector3(0, 0, 0);

  double phi = atan2(y, x);
  if (phi < 0)
    phi += utl::kTwoPi;

  const double phiRef = fDiskPhase1;
  int iBest = -2;
  double bestDist = -1;
  for (int i = -1; i <= 1; ++i) {
    const double pphi = phi - phiRef + i*utl::kTwoPi;
    const double rr = rRef*exp(pphi * fTanPitch);
    if (bestDist < 0 || std::abs(r-rr) < bestDist) {
      bestDist =  std::abs(r-rr);
      iBest = i;
    }
  }
  if (iBest == 0) {
    const double phi0 = phi - log(r/rRef) / fTanPitch;

    // Eq. (16)
    const double deltaPhi0 = utl::DeltaPhi(phiRef, phi0);
    const double delta = deltaPhi0 / fSpurWidth;
    const double B = fDiskB1 * exp(-0.5*pow(delta, 2));

    // Eq. (18)
    const double wS = 5*utl::degree;
    const double phiC = fSpurCenter;
    const double deltaPhiC = utl::DeltaPhi(phiC, phi);
    const double lC = fSpurLength;
    const double gS = 1 - utl::Sigmoid(std::abs(deltaPhiC), lC, wS);

    // Eq. (13)
    const double hd = 1 - utl::Sigmoid(std::abs(z), fDiskH, fDiskW);

    // Eq. (17)
    const double bS = rRef/r * B * hd * gS;
    const double bCyl[3] = {bS * fSinPitch, bS * fCosPitch, 0};
    const double cosPhi = x / r;
    const double sinPhi = y / r;
    return utl::CylToCart(bCyl, cosPhi, sinPhi);
  }
  else
    return Vector3(0, 0, 0);

}

vector UFMagneticField::GetSpiralField(const double x, const double y, const double z, const UFMagneticField &p)
  const
{
  // reference radius
  const double rRef = 5*utl::kpc;
  // inner boundary of spiral field
  const double rInner = 5*utl::kpc;
  const double wInner = 0.5*utl::kpc;
  // outer boundary of spiral field
  const double rOuter = 20*utl::kpc;
  const double wOuter = 0.5*utl::kpc;

  // cylindrical coordinates
  const double r2 = x*x + y*y;
  if (r2 == 0)
    return Vector3(0, 0, 0);
  const double r = sqrt(r2);
  const double phi = atan2(y, x);

  // Eq.(13)
  const double hdz = 1 - utl::Sigmoid(std::abs(z), fDiskH, fDiskW);

  // Eq.(14) times rRef divided by r
  const double rFacI = utl::Sigmoid(r, rInner, wInner);
  const double rFacO = 1 - utl::Sigmoid(r, rOuter, wOuter);
  // (using lim r--> 0 (1-exp(-r^2))/r --> r - r^3/2 + ...)
  const double rFac =  r > 1e-5*utl::pc ? (1-exp(-r*r)) / r : r * (1 - r2/2);
  const double gdrTimesRrefByR = rRef * rFac * rFacO * rFacI;

  // Eq. (12)
  const double phi0 = phi - log(r/rRef) / fTanPitch;

  // Eq. (10)
  const double b =
    fDiskB1 * cos(1 * (phi0 - fDiskPhase1)) +
    fDiskB2 * cos(2 * (phi0 - fDiskPhase2)) +
    fDiskB3 * cos(3 * (phi0 - fDiskPhase3));

  // Eq. (11)
  const double fac = hdz * gdrTimesRrefByR;
  const double bCyl[3] =
    { b * fac * fSinPitch,
      b * fac * fCosPitch,
      0};

  const double cosPhi = x / r;
  const double sinPhi = y / r;
  return utl::CylToCart(bCyl, cosPhi, sinPhi);
}
