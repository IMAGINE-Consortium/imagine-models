#ifndef REGULAR_TRAMPOLINE_H
#define REGULAR_TRAMPOLINE_H

#include "ImagineModels/RegularField.h"

class PyRegularVectorField : public RegularVectorField {
public:
    using RegularVectorField::RegularVectorField;
    vector at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(vector, RegularVectorField, at_position, x, y, z); }
};

class PyRegularScalarField : public RegularScalarField {
public:
    using RegularScalarField::RegularScalarField;
    number at_position(const double& x, const double& y, const double& z) const override {PYBIND11_OVERRIDE_PURE(number, RegularScalarField, at_position, x, y, z); }
};

#endif
