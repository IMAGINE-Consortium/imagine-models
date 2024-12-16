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

#include "UngerFarrar.h"
#include "helpers.h"
#include "units.h"


vector UFMagneticField::_at_position(const double &x, const double &y, const double &z, const UFMagneticField &p) const
{
  vector B_cart{{0., 0., 0.}};
  double squared_length = pow(x, 2) + pow(y, 2) + pow(z, 2);
  if (squared_length > p.fMaxRadiusSquared)
    return B_cart;
  else {
    const auto diskField = GetDiskField(x, y, z, p);
    const auto haloField = GetHaloField(x, y, z, p);
    return (diskField + haloField);
  }
}

vector UFMagneticField::GetDiskField(const double &x, const double &y, const double &z, const UFMagneticField &p)
  const
{
  if (p.activeModel == "spur")
    return GetSpurField(x, y, z, p);
  else
    return GetSpiralField(x, y, z, p);
}


vector UFMagneticField::GetHaloField(const double &x, const double &y, const double &z, const UFMagneticField &p)
  const
{
  if (p.activeModel == "twistX")
    return GetTwistedHaloField(x, y, z, p);
  else
    return
      GetToroidalHaloField(x, y, z, p) +
      GetPoloidalHaloField(x, y, z, p);
}


vector UFMagneticField::GetTwistedHaloField(const double x, const double y, const double z, const UFMagneticField &p)
  const
{
  const double r = sqrt(x*x + y*y);
  const double cosPhi = r > std::numeric_limits<double>::min() ? x / r : 1;
  const double sinPhi = r > std::numeric_limits<double>::min() ? y / r : 0;

  vector bXCart = GetPoloidalHaloField(x, y, z, p);
  vector bXCartTmp{{bXCart[0], bXCart[1], bXCart[2]}};
  vector bXCyl = Cart2Cyl(bXCartTmp, cosPhi, sinPhi);

  number bZ = bXCyl[2];
  number bR = bXCyl[0];

  number bPhi = 0;

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
  vector bCylX{{bR, bPhi , bZ}};
  return Cyl2Cart<vector>(bCylX, cosPhi, sinPhi);
}

vector UFMagneticField::GetToroidalHaloField(const double x, const double y, const double z, const UFMagneticField &p)
  const
{
  const double r2 = x*x + y*y;
  const double r = sqrt(r2);
  const double absZ = abs(z);

  number b0 = z >= 0 ? p.fToroidalBN : p.fToroidalBS;
  number rh = p.fToroidalR;
  number z0 = p.fToroidalZ;
  number fwh = p.fToroidalW;
  //number sigmoidR = Sigmoid<number>(r, rh, fwh);
  number sigmoidR = 1 / (1 + exp(-(r-rh)/fwh));
  //number sigmoidZ = Sigmoid<number>(absZ, p.fDiskH, p.fDiskW);
  number sigmoidZ = 1 / (1 + exp(-(absZ-p.fDiskH)/p.fDiskW));

  // Eq. (21)
  number bPhi = b0 * (1. - sigmoidR) * sigmoidZ * exp(-absZ/z0);

  vector bCyl{{0., bPhi, 0.}};
  const double cosPhi = r > std::numeric_limits<double>::min() ? x / r : 1;
  const double sinPhi = r > std::numeric_limits<double>::min() ? y / r : 0;
  return Cyl2Cart<vector>(bCyl, cosPhi, sinPhi);
}

vector UFMagneticField::GetPoloidalHaloField(const double x, const double y, const double z, const UFMagneticField &p)
  const
{
  const double r2 = x*x + y*y;
  const double r = std::sqrt(r2);

  number c = pow(p.fPoloidalA/p.fPoloidalZ, p.fPoloidalP);
  number a0p = pow(p.fPoloidalA, p.fPoloidalP);
  number rp = pow(r, p.fPoloidalP);
  number abszp = pow(abs(z), p.fPoloidalP);
  number cabszp = c*abszp;

  /*
    since $\sqrt{a^2 + b} - a$ is numerical unstable for $b\ll a$,
    we use $(\sqrt{a^2 + b} - a) \frac{\sqrt{a^2 + b} + a}{\sqrt{a^2
    + b} + a} = \frac{b}{\sqrt{a^2 + b} + a}$}
  */

  number t0 = a0p + cabszp - rp;
  number t1 = sqrt(pow(t0, 2) + 4*a0p*rp);
  number ap = 2*a0p*rp / (t1  + t0);

  number a = 0;
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
  number radialDependence =
    p.activeModel == "base" ?
    exp(-a/p.fPoloidalR) :
    //1 - Sigmoid<number>(a, p.fPoloidalR, p.fPoloidalW);
    1 - 1 / (1 + exp(-(a-p.fPoloidalR)/p.fPoloidalW));

  // Eq.(28)
  number Bzz = p.fPoloidalB * radialDependence;

  // (r/a)
  number rOverA =  1 / pow(2*a0p / (t1  + t0), 1/p.fPoloidalP);

  // Eq.(35) for p=n
  const double signZ = z < 0 ? -1 : 1;
  number Br =
    Bzz * c * a / rOverA * signZ * pow(abs(z), p.fPoloidalP - 1) / t1;

  // Eq.(36) for p=n
  number Bz = Bzz * pow(rOverA, p.fPoloidalP-2) * (ap + a0p) / t1;

  if (r < std::numeric_limits<double>::min())
    return vector{{0., 0., Bz}};
  else {
    vector bCylX{{Br, 0 , Bz}};
    const double cosPhi =  x / r;
    const double sinPhi =  y / r;
    return Cyl2Cart<vector>(bCylX, cosPhi, sinPhi);
  }
}

vector UFMagneticField::GetSpurField(const double x, const double y, const double z, const UFMagneticField &p)
  const
{
  // reference approximately at solar radius
  const double rRef = 8.2; //kpc

  // cylindrical coordinates
  const double r2 = x*x + y*y;
  const double r = sqrt(r2);
  if (r < std::numeric_limits<double>::min())
    return vector{{0, 0, 0}};

  double phi = atan2(y, x);
  if (phi < 0)
    phi += num::twopi;

  number phiRef = p.fDiskPhase1;
  int iBest = -2;
  number bestDist = -1;
  for (int i = -1; i <= 1; ++i) {
    number pphi = phi - phiRef + i*num::twopi;
    number rr = rRef*exp(pphi * p.fTanPitch);
    if (bestDist < 0 || abs(r-rr) < bestDist) {
      bestDist =  abs(r-rr);
      iBest = i;
    }
  }
  if (iBest == 0) {
    number phi0 = phi - log(r/rRef) / p.fTanPitch;

    // Eq. (16)
    //number deltaPhi0 = DeltaPhi<number, number, number>(phiRef, phi0);
    number deltaPhi0 = acos(cos(phi0)*cos(phiRef) + sin(phi0)*sin(phiRef));
    number delta = deltaPhi0 / p.fSpurWidth;
    number B = p.fDiskB1 * exp(-0.5*pow(delta, 2));

    // Eq. (18)
    const double wS = 5*num::rad;
    number phiC = p.fSpurCenter;
    //number deltaPhiC = DeltaPhi<number, number, number>(phiC, phi);
    number deltaPhiC = acos(cos(phi)*cos(phiC) + sin(phi)*sin(phiC));
    number lC = p.fSpurLength;
    //number gS = 1 - Sigmoid<number>(abs(deltaPhiC), lC, wS);
    number gS = 1 - 1 / (1 + exp(-(abs(deltaPhiC)-lC)/wS));

    // Eq. (13)
    //number hd = 1 - Sigmoid<number>(abs(z), p.fDiskH, p.fDiskW);
    number hd = 1 - 1 / (1 + exp(-(abs(z)-p.fDiskH)/p.fDiskW));

    // Eq. (17)
    number bS = rRef/r * B * hd * gS;
    vector bCyl{{bS * fSinPitch, bS * fCosPitch, 0.}};
    const double cosPhi = x / r;
    const double sinPhi = y / r;
    return Cyl2Cart<vector>(bCyl, cosPhi, sinPhi);
  }
  else
    return vector{{0, 0, 0}};

}

vector UFMagneticField::GetSpiralField(const double x, const double y, const double z, const UFMagneticField &p)
  const
{
  // reference radius
  const double rRef = 5.; // kpc
  // inner boundary of spiral field
  const double rInner = 5; // kpc
  const double wInner = 0.5; // kpc
  // outer boundary of spiral field
  const double rOuter = 20; // kpc
  const double wOuter = 0.5; // kpc

  // cylindrical coordinates
  const double r2 = x*x + y*y;
  if (r2 == 0)
    return vector{{0, 0, 0}};

  const double r = std::sqrt(r2);
  const double phi = std::atan2(y, x);

  // Eq.(13)
  //number hdz = 1 - Sigmoid(abs(z), fDiskH, fDiskW);
  number hdz = 1 - 1 / (1 + exp(-(abs(z)-p.fDiskH)/p.fDiskW));

  // Eq.(14) times rRef divided by r
  //const double rFacI = Sigmoid(r, rInner, wInner);
  const double rFacI = 1 / (1 + exp(-(r-rInner)/wInner));
  //const double rFacO = 1 - Sigmoid(r, rOuter, wOuter);
  const double rFacO = 1 - 1 / (1 + exp(-(r-rOuter)/wOuter));
  
  // (using lim r--> 0 (1-exp(-r^2))/r --> r - r^3/2 + ...)
  const double rFac =  r > 1e-5*astro::pc ? (1-exp(-r*r)) / r : r * (1 - r2/2);
  const double gdrTimesRrefByR = rRef * rFac * rFacO * rFacI;

  // Eq. (12)
  number phi0 = phi - log(r/rRef) / p.fTanPitch;

  // Eq. (10)
  number b =
    fDiskB1 * cos(1 * (phi0 - p.fDiskPhase1)) +
    fDiskB2 * cos(2 * (phi0 - p.fDiskPhase2)) +
    fDiskB3 * cos(3 * (phi0 - p.fDiskPhase3));

  // Eq. (11)
  number fac = hdz * gdrTimesRrefByR;
  vector bCyl{{ b * fac * p.fSinPitch,
      b * fac * p.fCosPitch,
      0.}};

  const double cosPhi = x / r;
  const double sinPhi = y / r;
  return Cyl2Cart<vector>(bCyl, cosPhi, sinPhi);
}
