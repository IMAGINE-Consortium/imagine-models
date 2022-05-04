#include "../c_library/headers/MagneticField.h"

// These classes are necessary to override virtual functions when binding

class PyMagneticField : public MagneticField {
public:
    using MagneticField::MagneticField; // Inherit constructors
};
class PyRegularMagneticField : public RegularMagneticField {
public:
    using RegularMagneticField::RegularMagneticField; // Inherit constructors
    std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE_PURE(std::vector<double>, RegularMagneticField, evaluate_at_pos, pos); }
};

class PyJF12MagneticField : public JF12MagneticField {
public:
    using JF12MagneticField::JF12MagneticField; // Inherit constructors
    std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE(std::vector<double>, JF12MagneticField, evaluate_at_pos, pos); }
};

class PyHelixMagneticField : public HelixMagneticField {
public:
    using HelixMagneticField::HelixMagneticField; // Inherit constructors
    std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE(std::vector<double>, HelixMagneticField, evaluate_at_pos, pos); }

};

class PyJaffeMagneticField : public JaffeMagneticField {
public:
    using JaffeMagneticField::JaffeMagneticField; // Inherit constructors
    std::vector<double> evaluate_at_pos(const std::vector<double> &pos) const override {PYBIND11_OVERRIDE(std::vector<double>, JaffeMagneticField, evaluate_at_pos, pos); }

};
