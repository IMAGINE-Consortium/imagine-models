#ifndef REGULAR_TRAMPOLINE_H
#define REGULAR_TRAMPOLINE_H

#include "ImagineModels/RegularField.h"

class PyRegularVectorField : public RegularVectorField {
public:
    using RegularVectorField::RegularVectorField;
    vector at_position(const double& x, const double& y, const double& z) const override {
        py::gil_scoped_acquire gil;
        py::function override = py::get_override(static_cast<const RegularVectorField *>(this), "at_position");
        if (!override)
            py::pybind11_fail("Tried to call pure virtual function \"RegularVectorField::at_position\"");
        auto b = override(x, y, z).cast<std::array<double, 3>>();
        return vector{{b[0], b[1], b[2]}};
    }
};

class PyRegularScalarField : public RegularScalarField {
public:
    using RegularScalarField::RegularScalarField;
    number at_position(const double& x, const double& y, const double& z) const override {
        py::gil_scoped_acquire gil;
        py::function override = py::get_override(static_cast<const RegularScalarField *>(this), "at_position");
        if (!override)
            py::pybind11_fail("Tried to call pure virtual function \"RegularScalarField::at_position\"");
        return number(override(x, y, z).cast<double>());
    }
};

#endif
