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
#include <map>
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

  double fMaxRadius = 20;

  number fPoloidalA    =  1 * astro::gpc;

  /// model parameters, see Table 3 of UF23 paper

  number fDiskB1        =  1.0878565e+00 * astro::microgauss;
  number fDiskB2        =  2.6605034e+00 * astro::microgauss;
  number fDiskB3        =  3.1166311e+00 * astro::microgauss;
  number fDiskH         =  7.9408965e-01 * astro::kpc;
  number fDiskPhase1    =  2.6316589e+02 * num::rad;
  number fDiskPhase2    =  9.7782269e+01 * num::rad;
  number fDiskPhase3    =  3.5112281e+01 * num::rad;
  number fDiskPitch     =  1.0106900e+01 * num::rad;
  number fDiskW         =  1.0720909e-01 * astro::kpc;
  number fPoloidalB     =  9.7775487e-01 * astro::microgauss;
  number fPoloidalP     =  1.4266186e+00 * astro::kpc;
  number fPoloidalR     =  7.2925417e+00 * astro::kpc;
  number fPoloidalW     =  1.1188158e-01 * astro::kpc;
  number fPoloidalZ     =  4.4597373e+00 * astro::kpc;
  number fStriation     =  3.4557571e-01;
  number fToroidalBN    =  3.2556760e+00 * astro::microgauss;
  number fToroidalBS    = -3.0914569e+00 * astro::microgauss;
  number fToroidalR     =  1.0193815e+01 * astro::kpc;
  number fToroidalW     =  1.6936993e+00 * astro::kpc;
  number fToroidalZ     =  4.0242749e+00 * astro::kpc;
  
  number fSpurCenter = 0;
  number fSpurLength = 0;
  number fSpurWidth = 0;
  number fTwistingTime = 0;

  std::map<std::string, std::map<std::string, double>> all_parameters = 
    {{"base", {
      {"fDiskB1"        ,  1.0878565e+00 * astro::microgauss},
      {"fDiskB2"        ,  2.6605034e+00 * astro::microgauss},
      {"fDiskB3"        ,  3.1166311e+00 * astro::microgauss},
      {"fDiskH"         ,  7.9408965e-01 * astro::kpc},
      {"fDiskPhase1"    ,  2.6316589e+02 * num::rad},
      {"fDiskPhase2"    ,  9.7782269e+01 * num::rad},
      {"fDiskPhase3"    ,  3.5112281e+01 * num::rad},
      {"fDiskPitch"     ,  1.0106900e+01 * num::rad},
      {"fDiskW"         ,  1.0720909e-01 * astro::kpc},
      {"fPoloidalB"     ,  9.7775487e-01 * astro::microgauss},
      {"fPoloidalP"     ,  1.4266186e+00 * astro::kpc},
      {"fPoloidalR"     ,  7.2925417e+00 * astro::kpc},
      {"fPoloidalW"     ,  1.1188158e-01 * astro::kpc},
      {"fPoloidalZ"     ,  4.4597373e+00 * astro::kpc},
      {"fStriation"     ,  3.4557571e-01},
      {"fToroidalBN"    ,  3.2556760e+00 * astro::microgauss},
      {"fToroidalBS"    , -3.0914569e+00 * astro::microgauss},
      {"fToroidalR"     ,  1.0193815e+01 * astro::kpc},
      {"fToroidalW"     ,  1.6936993e+00 * astro::kpc},
      {"fToroidalZ"     ,  4.0242749e+00 * astro::kpc},
      }}, 
     {"neCL", {
      {"fDiskB1"        ,  1.4259645e+00 * astro::microgauss},
      {"fDiskB2"        ,  1.3543223e+00 * astro::microgauss},
      {"fDiskB3"        ,  3.4390669e+00 * astro::microgauss},
      {"fDiskH"         ,  6.7405199e-01 * astro::kpc},
      {"fDiskPhase1"    ,  1.9961898e+02 * num::rad},
      {"fDiskPhase2"    ,  1.3541461e+02 * num::rad},
      {"fDiskPhase3"    ,  6.4909767e+01 * num::rad},
      {"fDiskPitch"     ,  1.1867859e+01 * num::rad},
      {"fDiskW"         ,  6.1162799e-02 * astro::kpc},
      {"fPoloidalB"     ,  9.8387831e-01 * astro::microgauss},
      {"fPoloidalP"     ,  1.6773615e+00 * astro::kpc},
      {"fPoloidalR"     ,  7.4084361e+00 * astro::kpc},
      {"fPoloidalW"     ,  1.4168192e-01 * astro::kpc},
      {"fPoloidalZ"     ,  3.6521188e+00 * astro::kpc},
      {"fStriation"     ,  3.3600213e-01},
      {"fToroidalBN"    ,  2.6256593e+00 * astro::microgauss},
      {"fToroidalBS"    , -2.5699466e+00 * astro::microgauss},
      {"fToroidalR"     ,  1.0134257e+01 * astro::kpc},
      {"fToroidalW"     ,  1.1547728e+00 * astro::kpc},
      {"fToroidalZ"     ,  4.5585463e+00 * astro::kpc},
     }}, 
     {"expX", {
    {"fDiskB1"        ,  9.9258148e-01 * astro::microgauss},
    {"fDiskB2"        ,  2.1821124e+00 * astro::microgauss},
    {"fDiskB3"        ,  3.1197345e+00 * astro::microgauss},
    {"fDiskH"         ,  7.1508681e-01 * astro::kpc},
    {"fDiskPhase1"    ,  2.4745741e+02 * num::rad},
    {"fDiskPhase2"    ,  9.8578879e+01 * num::rad},
    {"fDiskPhase3"    ,  3.4884485e+01 * num::rad},
    {"fDiskPitch"     ,  1.0027070e+01 * num::rad},
    {"fDiskW"         ,  9.8524736e-02 * astro::kpc},
    {"fPoloidalA"     ,  6.1938701e+00 * astro::kpc},
    {"fPoloidalB"     ,  5.8357990e+00 * astro::microgauss},
    {"fPoloidalP"     ,  1.9510779e+00 * astro::kpc},
    {"fPoloidalR"     ,  2.4994376e+00 * astro::kpc},
    {"fPoloidalZ"     ,  2.3684453e+00 * astro::kpc},
    {"fStriation"     ,  5.1440500e-01},
    {"fToroidalBN"    ,  2.7077434e+00 * astro::microgauss},
    {"fToroidalBS"    , -2.5677104e+00 * astro::microgauss},
    {"fToroidalR"     ,  1.0134022e+01 * astro::kpc},
    {"fToroidalW"     ,  2.0956159e+00 * astro::kpc},
    {"fToroidalZ"     ,  5.4564991e+00 * astro::kpc},
     }}, 
     {"spur", {
    {"fDiskB1"        , -4.2993328e+00 * astro::microgauss},
    {"fDiskH"         ,  7.5019749e-01 * astro::kpc},
    {"fDiskPhase1"    ,  1.5589875e+02 * num::rad},
    {"fDiskPitch"     ,  1.2074432e+01 * num::rad},
    {"fDiskW"         ,  1.2263120e-01 * astro::kpc},
    {"fPoloidalB"     ,  9.9302987e-01 * astro::microgauss},
    {"fPoloidalP"     ,  1.3982374e+00 * astro::kpc},
    {"fPoloidalR"     ,  7.1973387e+00 * astro::kpc},
    {"fPoloidalW"     ,  1.2262244e-01 * astro::kpc},
    {"fPoloidalZ"     ,  4.4853270e+00 * astro::kpc},
    {"fSpurCenter"    ,  1.5718686e+02 * num::rad},
    {"fSpurLength"    ,  3.1839577e+01 * num::rad},
    {"fSpurWidth"     ,  1.0318114e+01 * num::rad},
    {"fStriation"     ,  3.3022369e-01},
    {"fToroidalBN"    ,  2.9286724e+00 * astro::microgauss},
    {"fToroidalBS"    , -2.5979895e+00 * astro::microgauss},
    {"fToroidalR"     ,  9.7536425e+00 * astro::kpc},
    {"fToroidalW"     ,  1.4210055e+00 * astro::kpc},
    {"fToroidalZ"     ,  6.0941229e+00 * astro::kpc},
     }}, 
     {"cre10", {
    {"fDiskB1"        ,  1.2035697e+00 * astro::microgauss},
    {"fDiskB2"        ,  2.7478490e+00 * astro::microgauss},
    {"fDiskB3"        ,  3.2104342e+00 * astro::microgauss},
    {"fDiskH"         ,  8.0844932e-01 * astro::kpc},
    {"fDiskPhase1"    ,  2.6515882e+02 * num::rad},
    {"fDiskPhase2"    ,  9.8211313e+01 * num::rad},
    {"fDiskPhase3"    ,  3.5944588e+01 * num::rad},
    {"fDiskPitch"     ,  1.0162759e+01 * num::rad},
    {"fDiskW"         ,  1.0824003e-01 * astro::kpc},
    {"fPoloidalB"     ,  9.6938453e-01 * astro::microgauss},
    {"fPoloidalP"     ,  1.4150957e+00 * astro::kpc},
    {"fPoloidalR"     ,  7.2987296e+00 * astro::kpc},
    {"fPoloidalW"     ,  1.0923051e-01 * astro::kpc},
    {"fPoloidalZ"     ,  4.5748332e+00 * astro::kpc},
    {"fStriation"     ,  2.4950386e-01},
    {"fToroidalBN"    ,  3.7308133e+00 * astro::microgauss},
    {"fToroidalBS"    , -3.5039958e+00 * astro::microgauss},
    {"fToroidalR"     ,  1.0407507e+01 * astro::kpc},
    {"fToroidalW"     ,  1.7398375e+00 * astro::kpc},
    {"fToroidalZ"     ,  2.9272800e+00 * astro::kpc},
     }}, 
     {"synCG", {
    {"fDiskB1"        ,  8.1386878e-01 * astro::microgauss},
    {"fDiskB2"        ,  2.0586930e+00 * astro::microgauss},
    {"fDiskB3"        ,  2.9437335e+00 * astro::microgauss},
    {"fDiskH"         ,  6.2172353e-01 * astro::kpc},
    {"fDiskPhase1"    ,  2.2988551e+02 * num::rad},
    {"fDiskPhase2"    ,  9.7388282e+01 * num::rad},
    {"fDiskPhase3"    ,  3.2927367e+01 * num::rad},
    {"fDiskPitch"     ,  9.9034844e+00 * num::rad},
    {"fDiskW"         ,  6.6517521e-02 * astro::kpc},
    {"fPoloidalB"     ,  8.0883734e-01 * astro::microgauss},
    {"fPoloidalP"     ,  1.5820957e+00 * astro::kpc},
    {"fPoloidalR"     ,  7.4625235e+00 * astro::kpc},
    {"fPoloidalW"     ,  1.5003765e-01 * astro::kpc},
    {"fPoloidalZ"     ,  3.5338550e+00 * astro::kpc},
    {"fStriation"     ,  6.3434763e-01},
    {"fToroidalBN"    ,  2.3991193e+00 * astro::microgauss},
    {"fToroidalBS"    , -2.0919944e+00 * astro::microgauss},
    {"fToroidalR"     ,  9.4227834e+00 * astro::kpc},
    {"fToroidalW"     ,  9.1608418e-01 * astro::kpc},
    {"fToroidalZ"     ,  5.5844594e+00 * astro::kpc},
     }}, 
     {"twistX", {
    {"fDiskB1"        ,  1.3741995e+00 * astro::microgauss},
    {"fDiskB2"        ,  2.0089881e+00 * astro::microgauss},
    {"fDiskB3"        ,  1.5212463e+00 * astro::microgauss},
    {"fDiskH"         ,  9.3806180e-01 * astro::kpc},
    {"fDiskPhase1"    ,  2.3560316e+02 * num::rad},
    {"fDiskPhase2"    ,  1.0189856e+02 * num::rad},
    {"fDiskPhase3"    ,  5.6187572e+01 * num::rad},
    {"fDiskPitch"     ,  1.2100979e+01 * num::rad},
    {"fDiskW"         ,  1.4933338e-01 * astro::kpc},
    {"fPoloidalB"     ,  6.2793114e-01 * astro::microgauss},
    {"fPoloidalP"     ,  2.3292519e+00 * astro::kpc},
    {"fPoloidalR"     ,  7.9212358e+00 * astro::kpc},
    {"fPoloidalW"     ,  2.9056201e-01 * astro::kpc},
    {"fPoloidalZ"     ,  2.6274437e+00 * astro::kpc},
    {"fStriation"     ,  7.7616317e-01},
    {"fTwistingTime"  ,  5.4733549e+01 * astro::megayear},
     }}, 
     {"nebCor", {
    {"fDiskB1"        ,  1.4081935e+00 * astro::microgauss},
    {"fDiskB2"        ,  3.5292400e+00 * astro::microgauss},
    {"fDiskB3"        ,  4.1290147e+00 * astro::microgauss},
    {"fDiskH"         ,  8.1151971e-01 * astro::kpc},
    {"fDiskPhase1"    ,  2.6447529e+02 * num::rad},
    {"fDiskPhase2"    ,  9.7572660e+01 * num::rad},
    {"fDiskPhase3"    ,  3.6403798e+01 * num::rad},
    {"fDiskPitch"     ,  1.0151183e+01 * num::rad},
    {"fDiskW"         ,  1.1863734e-01 * astro::kpc},
    {"fPoloidalB"     ,  1.3485916e+00 * astro::microgauss},
    {"fPoloidalP"     ,  1.3414395e+00 * astro::kpc},
    {"fPoloidalR"     ,  7.2473841e+00 * astro::kpc},
    {"fPoloidalW"     ,  1.4318227e-01 * astro::kpc},
    {"fPoloidalZ"     ,  4.8242603e+00 * astro::kpc},
    {"fStriation"     ,  3.8610837e-10},
    {"fToroidalBN"    ,  4.6491142e+00 * astro::microgauss},
    {"fToroidalBS"    , -4.5006610e+00 * astro::microgauss},
    {"fToroidalR"     ,  1.0205288e+01 * astro::kpc},
    {"fToroidalW"     ,  1.7004868e+00 * astro::kpc},
    {"fToroidalZ"     ,  3.5557767e+00 * astro::kpc},
     }}, 
    };

  vector at_position(const double &x, const double &y, const double &z) const
  {
    return _at_position(x, y, z, *this);
  }

  void set_parameters(const std::string &model_choice);


#if autodiff_FOUND
  const std::set<std::string> all_diff{"fDiskB1", "fDiskB2", "fDiskB3", "fDiskH", "fDiskPhase1", "fDiskPhase2", "fDiskPhase3", "fDiskPitch", "fDiskW", "fPoloidalA", "fPoloidalB", "fPoloidalP", "fPoloidalR", "fPoloidalW", "fPoloidalZ", "fSpurCenter", "fSpurLength", "fSpurWidth", "fStriation", "fToroidalBN", "fToroidalBS", "fToroidalR", "fToroidalW", "fToroidalZ", "fTwistingTime"};
  std::set<std::string> active_diff{"fDiskB1", "fDiskB2", "fDiskB3", "fDiskH", "fDiskPhase1", "fDiskPhase2", "fDiskPhase3", "fDiskPitch", "fDiskW", "fPoloidalA", "fPoloidalB", "fPoloidalP", "fPoloidalR", "fPoloidalW", "fPoloidalZ", "fSpurCenter", "fSpurLength", "fSpurWidth", "fStriation", "fToroidalBN", "fToroidalBS", "fToroidalR", "fToroidalW", "fToroidalZ", "fTwistingTime"};

  Eigen::MatrixXd derivative(const double &x, const double &y, const double &z)
  {
    return _jac(x, y, z, *this);
  }
#endif

private:

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
