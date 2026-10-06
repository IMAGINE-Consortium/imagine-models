#ifndef RANDOM_TRAMPOLINE_H
#define RANDOM_TRAMPOLINE_H

#include "ImagineModelsRandom/RandomScalarField.h"
#include "ImagineModelsRandom/RandomVectorField.h"

class PyRandomVectorField : public RandomVectorField {
public:
    using RandomVectorField::RandomVectorField;
    double spatial_profile(const double &x, const double &y, const double &z) const override {PYBIND11_OVERRIDE_PURE(double, RandomVectorField, spatial_profile, x, y, z); }
    double calculate_fourier_sigma(const double &abs_k, const double &dk) const override {PYBIND11_OVERRIDE_PURE(double, RandomVectorField, calculate_fourier_sigma, abs_k, dk); }
};

class PyRandomScalarField : public RandomScalarField {
public:
    using RandomScalarField::RandomScalarField;
    double spatial_profile(const double &x, const double &y, const double &z) const override {PYBIND11_OVERRIDE_PURE(double, RandomScalarField, spatial_profile, x, y, z); }
    double calculate_fourier_sigma(const double &abs_k, const double &dk) const override {PYBIND11_OVERRIDE_PURE(double, RandomScalarField, calculate_fourier_sigma, abs_k, dk); }
};

#endif
