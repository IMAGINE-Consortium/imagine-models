#ifndef RANDOM_TRAMPOLINE_H
#define RANDOM_TRAMPOLINE_H

#include "ImagineModelsRandom/RandomScalarField.h"
#include "ImagineModelsRandom/RandomVectorField.h"

class PyRandomVectorField : public RandomVectorField {
public:
    using RandomVectorField::RandomVectorField;
    double spectrum(const double &abs_k) const override {PYBIND11_OVERRIDE_PURE(double, RandomVectorField, spectrum, abs_k); }
    double rms(const double &x, const double &y, const double &z) const override {PYBIND11_OVERRIDE_PURE(double, RandomVectorField, rms, x, y, z); }
    Vec3<double> anisotropy_direction(const double &x, const double &y, const double &z) const override {PYBIND11_OVERRIDE(Vec3<double>, RandomVectorField, anisotropy_direction, x, y, z); }
};

class PyRandomScalarField : public RandomScalarField {
public:
    using RandomScalarField::RandomScalarField;
    double spectrum(const double &abs_k) const override {PYBIND11_OVERRIDE_PURE(double, RandomScalarField, spectrum, abs_k); }
    double rms(const double &x, const double &y, const double &z) const override {PYBIND11_OVERRIDE_PURE(double, RandomScalarField, rms, x, y, z); }
    double mean(const double &x, const double &y, const double &z) const override {PYBIND11_OVERRIDE(double, RandomScalarField, mean, x, y, z); }
};

#endif
