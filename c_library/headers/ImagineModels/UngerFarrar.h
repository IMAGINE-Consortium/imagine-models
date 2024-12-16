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


#ifndef UNGERFARRAR_H
#define UNGERFARRAR_H

#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "RegularField.h"
#include "units.h"

class UFMagneticField : public RegularVectorField
{
protected:
  vector _at_position(const double &x, const double &y, const double &z, const UFMagneticField &p) const;

#if autodiff_FOUND
  Eigen::MatrixXd _jac(const double &x, const double &y, const double &z, UFMagneticField &p) const;
#endif

public:
  using RegularVectorField ::RegularVectorField;

public:
  /// model variations (see Tab.2 of UF23 paper)
  const std::array<std::string, 8> possibleModels{"base", "neCL", "expX", "spur", "cre10", "synCG", "twistX", "nebCor"};

  /// model type given in constructor
  std::string activeModel = "base";
  /// maximum galacto-centric radius beyond which B=0
  /// model parameters, see Table 3 of UF23 paper
  number fDiskB1 = 0;
  number fDiskB2 = 0;
  number fDiskB3 = 0;
  number fDiskH = 0;
  number fDiskPhase1 = 0;
  number fDiskPhase2 = 0;
  number fDiskPhase3 = 0;
  number fDiskPitch = 0;
  number fDiskW = 0;
  number fPoloidalA = 0;
  number fPoloidalB = 0;
  number fPoloidalP = 0;
  number fPoloidalR = 0;
  number fPoloidalW = 0;
  number fPoloidalZ = 0;
  number fSpurCenter = 0;
  number fSpurLength = 0;
  number fSpurWidth = 0;
  number fStriation = 0;
  number fToroidalBN = 0;
  number fToroidalBS = 0;
  number fToroidalR = 0;
  number fToroidalW = 0;
  number fToroidalZ = 0;
  number fTwistingTime = 0;


  vector at_position(const double &x, const double &y, const double &z) const
  {
    return _at_position(x, y, z, *this);
  }

#if autodiff_FOUND
  const std::set<std::string> all_diff{"fDiskB1", "fDiskB2", "fDiskB3", "fDiskH", "fDiskPhase1", "fDiskPhase2", "fDiskPhase3", "fDiskPitch", "fDiskW", "fPoloidalA", "fPoloidalB", "fPoloidalP", "fPoloidalR", "fPoloidalW", "fPoloidalZ", "fSpurCenter", "fSpurLength", "fSpurWidth", "fStriation", "fToroidalBN", "fToroidalBS", "fToroidalR", "fToroidalW", "fToroidalZ", "fTwistingTime"};
  std::set<std::string> active_diff{"fDiskB1", "fDiskB2", "fDiskB3", "fDiskH", "fDiskPhase1", "fDiskPhase2", "fDiskPhase3", "fDiskPitch", "fDiskW", "fPoloidalA", "fPoloidalB", "fPoloidalP", "fPoloidalR", "fPoloidalW", "fPoloidalZ", "fSpurCenter", "fSpurLength", "fSpurWidth", "fStriation", "fToroidalBN", "fToroidalBS", "fToroidalR", "fToroidalW", "fToroidalZ", "fTwistingTime"};

  Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
  {
    return _jac(x, y, z, *this);
  }
#endif

private:

  // some pre-calculated derived parameter values

  double fSinPitch = 0;
  double fCosPitch = 0;
  double fTanPitch = 0;
  const double fMaxRadiusSquared = 0;

  /// major field components
  vector GetDiskField(const double &x, const double &y, const double &z, const UFMagneticField &p) const;
  vector GetHaloField(const double &x, const double &y, const double &z, const UFMagneticField &p) const;

  /// sub-components depending on model type
  /// -- Sec. 5.2.2
  vector GetSpiralField(const double x, const double y, const double z, const UFMagneticField &p) const;
  /// -- Sec. 5.2.3
  vector GetSpurField(const double x, const double y, const double z, const UFMagneticField &p) const;
  /// -- Sec. 5.3.1
  vector GetToroidalHaloField(const double x, const double y, const double z, const UFMagneticField &p) const;
  /// -- Sec. 5.3.2
  vector GetPoloidalHaloField(const double x, const double y, const double z, const UFMagneticField &p) const;
  /// -- Sec. 5.3.3
  vector GetTwistedHaloField(const double x, const double y, const double z, const UFMagneticField &p) const;
  
};

#endif
