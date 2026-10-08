#ifndef REGULAR_TRAMPOLINE_H
#define REGULAR_TRAMPOLINE_H

#include "ImagineModels/RegularField.h"

class PyRegularVectorField : public RegularVectorField {
public:
    using RegularVectorField::RegularVectorField;
    Vec3<double> at_position(const double &x, const double &y, const double &z) const override {
        py::gil_scoped_acquire gil;
        py::function override = py::get_override(static_cast<const RegularVectorField *>(this), "at_position");
        if (!override)
            py::pybind11_fail("Tried to call pure virtual function \"RegularVectorField::at_position\"");
        auto b = override(x, y, z).cast<Vec3<double>>();
        return b;
    }
};

class PyRegularScalarField : public RegularScalarField {
public:
    using RegularScalarField::RegularScalarField;
    double at_position(const double &x, const double &y, const double &z) const override {
        py::gil_scoped_acquire gil;
        py::function override = py::get_override(static_cast<const RegularScalarField *>(this), "at_position");
        if (!override)
            py::pybind11_fail("Tried to call pure virtual function \"RegularScalarField::at_position\"");
        return override(x, y, z).cast<double>();
    }
};

#endif
