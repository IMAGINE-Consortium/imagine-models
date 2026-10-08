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

// Reference: Unger & Farrar 2024, arXiv:2311.12120
// Based on: authors' code (UF23Field v1.1, doi:10.5281/zenodo.11321212, BSD-2)

#pragma once

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>

#include "ImagineModels/RegularModel.h"
#include "ImagineModels/units.h"

namespace imagine {

#define UF24_PARAMETERS(X)                             \
    X(fPoloidalA, 1 * units::Gpc)                      \
    X(fDiskB1, 1.0878565e+00 * units::microgauss)      \
    X(fDiskB2, 2.6605034e+00 * units::microgauss)      \
    X(fDiskB3, 3.1166311e+00 * units::microgauss)      \
    X(fDiskH, 7.9408965e-01 * units::kpc)              \
    X(fDiskPhase1, 2.6316589e+02) /* deg */            \
    X(fDiskPhase2, 9.7782269e+01) /* deg */            \
    X(fDiskPhase3, 3.5112281e+01) /* deg */            \
    X(fDiskPitch, 1.0106900e+01)  /* deg */            \
    X(fDiskW, 1.0720909e-01 * units::kpc)              \
    X(fPoloidalB, 9.7775487e-01 * units::microgauss)   \
    X(fPoloidalP, 1.4266186e+00 * units::kpc)          \
    X(fPoloidalR, 7.2925417e+00 * units::kpc)          \
    X(fPoloidalW, 1.1188158e-01 * units::kpc)          \
    X(fPoloidalZ, 4.4597373e+00 * units::kpc)          \
    X(fStriation, 3.4557571e-01)                       \
    X(fToroidalBN, 3.2556760e+00 * units::microgauss)  \
    X(fToroidalBS, -3.0914569e+00 * units::microgauss) \
    X(fToroidalR, 1.0193815e+01 * units::kpc)          \
    X(fToroidalW, 1.6936993e+00 * units::kpc)          \
    X(fToroidalZ, 4.0242749e+00 * units::kpc)          \
    X(fSpurCenter, 0) /* deg */                        \
    X(fSpurLength, 0) /* deg */                        \
    X(fSpurWidth, 0)  /* deg */                        \
    X(fTwistingTime, 0)

IMAGINE_PARAMETERS(UF24Parameters, UF24_PARAMETERS)

class UF24MagneticField : public RegularVectorModel<UF24MagneticField, UF24Parameters> {
public:
    // variants, Table 2
    const std::array<std::string, 8> available_models{"base",  "neCL",  "expX",   "spur",
                                                      "cre10", "synCG", "twistX", "nebCor"};

    // B = 0 beyond this radius

    double fMaxRadius = 30;

    // parameters, Table 3

    std::map<std::string, std::map<std::string, double>> all_parameters = {
        {"base",
         {
             {"fDiskB1", 1.0878565e+00 * units::microgauss},
             {"fDiskB2", 2.6605034e+00 * units::microgauss},
             {"fDiskB3", 3.1166311e+00 * units::microgauss},
             {"fDiskH", 7.9408965e-01 * units::kpc},
             {"fDiskPhase1", 2.6316589e+02},
             {"fDiskPhase2", 9.7782269e+01},
             {"fDiskPhase3", 3.5112281e+01},
             {"fDiskPitch", 1.0106900e+01},
             {"fDiskW", 1.0720909e-01 * units::kpc},
             {"fPoloidalB", 9.7775487e-01 * units::microgauss},
             {"fPoloidalP", 1.4266186e+00 * units::kpc},
             {"fPoloidalR", 7.2925417e+00 * units::kpc},
             {"fPoloidalW", 1.1188158e-01 * units::kpc},
             {"fPoloidalZ", 4.4597373e+00 * units::kpc},
             {"fStriation", 3.4557571e-01},
             {"fToroidalBN", 3.2556760e+00 * units::microgauss},
             {"fToroidalBS", -3.0914569e+00 * units::microgauss},
             {"fToroidalR", 1.0193815e+01 * units::kpc},
             {"fToroidalW", 1.6936993e+00 * units::kpc},
             {"fToroidalZ", 4.0242749e+00 * units::kpc},
         }},
        {"neCL",
         {
             {"fDiskB1", 1.4259645e+00 * units::microgauss},
             {"fDiskB2", 1.3543223e+00 * units::microgauss},
             {"fDiskB3", 3.4390669e+00 * units::microgauss},
             {"fDiskH", 6.7405199e-01 * units::kpc},
             {"fDiskPhase1", 1.9961898e+02},
             {"fDiskPhase2", 1.3541461e+02},
             {"fDiskPhase3", 6.4909767e+01},
             {"fDiskPitch", 1.1867859e+01},
             {"fDiskW", 6.1162799e-02 * units::kpc},
             {"fPoloidalB", 9.8387831e-01 * units::microgauss},
             {"fPoloidalP", 1.6773615e+00 * units::kpc},
             {"fPoloidalR", 7.4084361e+00 * units::kpc},
             {"fPoloidalW", 1.4168192e-01 * units::kpc},
             {"fPoloidalZ", 3.6521188e+00 * units::kpc},
             {"fStriation", 3.3600213e-01},
             {"fToroidalBN", 2.6256593e+00 * units::microgauss},
             {"fToroidalBS", -2.5699466e+00 * units::microgauss},
             {"fToroidalR", 1.0134257e+01 * units::kpc},
             {"fToroidalW", 1.1547728e+00 * units::kpc},
             {"fToroidalZ", 4.5585463e+00 * units::kpc},
         }},
        {"expX",
         {
             {"fDiskB1", 9.9258148e-01 * units::microgauss},
             {"fDiskB2", 2.1821124e+00 * units::microgauss},
             {"fDiskB3", 3.1197345e+00 * units::microgauss},
             {"fDiskH", 7.1508681e-01 * units::kpc},
             {"fDiskPhase1", 2.4745741e+02},
             {"fDiskPhase2", 9.8578879e+01},
             {"fDiskPhase3", 3.4884485e+01},
             {"fDiskPitch", 1.0027070e+01},
             {"fDiskW", 9.8524736e-02 * units::kpc},
             {"fPoloidalA", 6.1938701e+00 * units::kpc},
             {"fPoloidalB", 5.8357990e+00 * units::microgauss},
             {"fPoloidalP", 1.9510779e+00 * units::kpc},
             {"fPoloidalR", 2.4994376e+00 * units::kpc},
             {"fPoloidalZ", 6.1938701e+00 * std::tan(2.0926122e+01 * units::deg) * units::kpc},
             {"fStriation", 5.1440500e-01},
             {"fToroidalBN", 2.7077434e+00 * units::microgauss},
             {"fToroidalBS", -2.5677104e+00 * units::microgauss},
             {"fToroidalR", 1.0134022e+01 * units::kpc},
             {"fToroidalW", 2.0956159e+00 * units::kpc},
             {"fToroidalZ", 5.4564991e+00 * units::kpc},
         }},
        {"spur",
         {
             {"fDiskB1", -4.2993328e+00 * units::microgauss},
             {"fDiskH", 7.5019749e-01 * units::kpc},
             {"fDiskPhase1", 1.5589875e+02},
             {"fDiskPitch", 1.2074432e+01},
             {"fDiskW", 1.2263120e-01 * units::kpc},
             {"fPoloidalB", 9.9302987e-01 * units::microgauss},
             {"fPoloidalP", 1.3982374e+00 * units::kpc},
             {"fPoloidalR", 7.1973387e+00 * units::kpc},
             {"fPoloidalW", 1.2262244e-01 * units::kpc},
             {"fPoloidalZ", 4.4853270e+00 * units::kpc},
             {"fSpurCenter", 1.5718686e+02},
             {"fSpurLength", 3.1839577e+01},
             {"fSpurWidth", 1.0318114e+01},
             {"fStriation", 3.3022369e-01},
             {"fToroidalBN", 2.9286724e+00 * units::microgauss},
             {"fToroidalBS", -2.5979895e+00 * units::microgauss},
             {"fToroidalR", 9.7536425e+00 * units::kpc},
             {"fToroidalW", 1.4210055e+00 * units::kpc},
             {"fToroidalZ", 6.0941229e+00 * units::kpc},
         }},
        {"cre10",
         {
             {"fDiskB1", 1.2035697e+00 * units::microgauss},
             {"fDiskB2", 2.7478490e+00 * units::microgauss},
             {"fDiskB3", 3.2104342e+00 * units::microgauss},
             {"fDiskH", 8.0844932e-01 * units::kpc},
             {"fDiskPhase1", 2.6515882e+02},
             {"fDiskPhase2", 9.8211313e+01},
             {"fDiskPhase3", 3.5944588e+01},
             {"fDiskPitch", 1.0162759e+01},
             {"fDiskW", 1.0824003e-01 * units::kpc},
             {"fPoloidalB", 9.6938453e-01 * units::microgauss},
             {"fPoloidalP", 1.4150957e+00 * units::kpc},
             {"fPoloidalR", 7.2987296e+00 * units::kpc},
             {"fPoloidalW", 1.0923051e-01 * units::kpc},
             {"fPoloidalZ", 4.5748332e+00 * units::kpc},
             {"fStriation", 2.4950386e-01},
             {"fToroidalBN", 3.7308133e+00 * units::microgauss},
             {"fToroidalBS", -3.5039958e+00 * units::microgauss},
             {"fToroidalR", 1.0407507e+01 * units::kpc},
             {"fToroidalW", 1.7398375e+00 * units::kpc},
             {"fToroidalZ", 2.9272800e+00 * units::kpc},
         }},
        {"synCG",
         {
             {"fDiskB1", 8.1386878e-01 * units::microgauss},
             {"fDiskB2", 2.0586930e+00 * units::microgauss},
             {"fDiskB3", 2.9437335e+00 * units::microgauss},
             {"fDiskH", 6.2172353e-01 * units::kpc},
             {"fDiskPhase1", 2.2988551e+02},
             {"fDiskPhase2", 9.7388282e+01},
             {"fDiskPhase3", 3.2927367e+01},
             {"fDiskPitch", 9.9034844e+00},
             {"fDiskW", 6.6517521e-02 * units::kpc},
             {"fPoloidalB", 8.0883734e-01 * units::microgauss},
             {"fPoloidalP", 1.5820957e+00 * units::kpc},
             {"fPoloidalR", 7.4625235e+00 * units::kpc},
             {"fPoloidalW", 1.5003765e-01 * units::kpc},
             {"fPoloidalZ", 3.5338550e+00 * units::kpc},
             {"fStriation", 6.3434763e-01},
             {"fToroidalBN", 2.3991193e+00 * units::microgauss},
             {"fToroidalBS", -2.0919944e+00 * units::microgauss},
             {"fToroidalR", 9.4227834e+00 * units::kpc},
             {"fToroidalW", 9.1608418e-01 * units::kpc},
             {"fToroidalZ", 5.5844594e+00 * units::kpc},
         }},
        {"twistX",
         {
             {"fDiskB1", 1.3741995e+00 * units::microgauss},
             {"fDiskB2", 2.0089881e+00 * units::microgauss},
             {"fDiskB3", 1.5212463e+00 * units::microgauss},
             {"fDiskH", 9.3806180e-01 * units::kpc},
             {"fDiskPhase1", 2.3560316e+02},
             {"fDiskPhase2", 1.0189856e+02},
             {"fDiskPhase3", 5.6187572e+01},
             {"fDiskPitch", 1.2100979e+01},
             {"fDiskW", 1.4933338e-01 * units::kpc},
             {"fPoloidalB", 6.2793114e-01 * units::microgauss},
             {"fPoloidalP", 2.3292519e+00 * units::kpc},
             {"fPoloidalR", 7.9212358e+00 * units::kpc},
             {"fPoloidalW", 2.9056201e-01 * units::kpc},
             {"fPoloidalZ", 2.6274437e+00 * units::kpc},
             {"fStriation", 7.7616317e-01},
             {"fTwistingTime", 5.4733549e+01 * units::megayear},
         }},
        {"nebCor",
         {
             {"fDiskB1", 1.4081935e+00 * units::microgauss},
             {"fDiskB2", 3.5292400e+00 * units::microgauss},
             {"fDiskB3", 4.1290147e+00 * units::microgauss},
             {"fDiskH", 8.1151971e-01 * units::kpc},
             {"fDiskPhase1", 2.6447529e+02},
             {"fDiskPhase2", 9.7572660e+01},
             {"fDiskPhase3", 3.6403798e+01},
             {"fDiskPitch", 1.0151183e+01},
             {"fDiskW", 1.1863734e-01 * units::kpc},
             {"fPoloidalB", 1.3485916e+00 * units::microgauss},
             {"fPoloidalP", 1.3414395e+00 * units::kpc},
             {"fPoloidalR", 7.2473841e+00 * units::kpc},
             {"fPoloidalW", 1.4318227e-01 * units::kpc},
             {"fPoloidalZ", 4.8242603e+00 * units::kpc},
             {"fStriation", 3.8610837e-10},
             {"fToroidalBN", 4.6491142e+00 * units::microgauss},
             {"fToroidalBS", -4.5006610e+00 * units::microgauss},
             {"fToroidalR", 1.0205288e+01 * units::kpc},
             {"fToroidalW", 1.7004868e+00 * units::kpc},
             {"fToroidalZ", 3.5557767e+00 * units::kpc},
         }},
    };

    explicit UF24MagneticField(const std::string &model = "base") { set_model(model); }

    void set_model(const std::string &model);
    const std::string &model() const { return active_model; }

private:
    // active variant
    std::string active_model = "base";

    // major field components
    template <typename T>
    Vec3<T> GetDiskField(const double &x, const double &y, const double &z, const UF24Parameters<T> &p) const;
    template <typename T>
    Vec3<T> GetHaloField(const double &x, const double &y, const double &z, const UF24Parameters<T> &p) const;

    // variant sub-components
    // -- Sec. 5.2.2
    template <typename T>
    Vec3<T> GetSpiralField(const double x, const double y, const double z, const UF24Parameters<T> &p) const;
    // -- Sec. 5.2.3
    template <typename T>
    Vec3<T> GetSpurField(const double x, const double y, const double z, const UF24Parameters<T> &p) const;
    // -- Sec. 5.3.1
    template <typename T>
    Vec3<T> GetToroidalHaloField(const double x, const double y, const double z, const UF24Parameters<T> &p) const;
    // -- Sec. 5.3.2
    template <typename T>
    Vec3<T> GetPoloidalHaloField(const double x, const double y, const double z, const UF24Parameters<T> &p) const;
    // -- Sec. 5.3.3
    template <typename T>
    Vec3<T> GetTwistedHaloField(const double x, const double y, const double z, const UF24Parameters<T> &p) const;

public:
    template <typename T>
    Vec3<T> field(const double &x, const double &y, const double &z, const UF24Parameters<T> &p) const;
};

}
