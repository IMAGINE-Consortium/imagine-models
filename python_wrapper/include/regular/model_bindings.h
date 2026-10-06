#ifndef MODEL_BINDINGS_H
#define MODEL_BINDINGS_H

#include <map>
#include <string>
#include <tuple>
#include <type_traits>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ImagineModels/RegularModel.h"

namespace py = pybind11;
using namespace pybind11::literals;

template <typename Model>
auto bind_regular_model(py::module_ &m, const char *name) {
    constexpr bool is_vector = std::is_base_of_v<RegularVectorField, Model>;
    using Base = std::conditional_t<is_vector, RegularVectorField, RegularScalarField>;
    py::class_<Model, Base> cls(m, name);
    cls.def(py::init<>());

    for (const std::string &parameter : Model::parameter_names())
        cls.def_property(parameter.c_str(),
            [parameter](const Model &self) { return self.get_parameter(parameter); },
            [parameter](Model &self, double value) { self.set_parameter(parameter, value); });

    cls.def_property("parameters", [](const Model &self) {
        py::dict out;
        for (const std::string &parameter : Model::parameter_names())
            out[parameter.c_str()] = self.get_parameter(parameter);
        return out;
    }, &Model::set_parameter_map);
    cls.def_property_readonly_static("parameter_names", [](py::object) { return Model::parameter_names(); });

#if IMAGINE_HAS_AUTODIFF
    cls.def_readwrite("active_parameters", &Model::active_parameters);
    cls.def("derivative", &Model::derivative, "x"_a, "y"_a, "z"_a);
#endif

    if constexpr (is_vector)
        cls.def("at_position", [](const Model &self, double x, double y, double z) {
            auto b = self.at_position(x, y, z);
            return std::make_tuple(static_cast<double>(b[0]), static_cast<double>(b[1]), static_cast<double>(b[2]));
        }, "x"_a, "y"_a, "z"_a);
    else
        cls.def("at_position", [](const Model &self, double x, double y, double z) {
            return static_cast<double>(self.at_position(x, y, z));
        }, "x"_a, "y"_a, "z"_a);

    return cls;
}

#endif
