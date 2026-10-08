#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "ImagineModels/Parameters.h"
#include "ImagineModels/RegularField.h"

namespace imagine {

template <typename Derived, template <typename> class Params> class RegularModel {
public:
    template <typename T> using parameters_t = Params<T>;

    Params<double> parameters;
    std::vector<std::string> active_parameters = parameter_names();

    static std::vector<std::string> parameter_names() {
        return std::vector<std::string>(Params<double>::names.begin(), Params<double>::names.end());
    }

    static std::size_t parameter_index(const std::string &name) {
        for (std::size_t i = 0; i < Params<double>::size; ++i)
            if (name == Params<double>::names[i])
                return i;
        throw std::invalid_argument("Unknown parameter '" + name + "'.");
    }

    double get_parameter(const std::string &name) const { return *parameters.addresses()[parameter_index(name)]; }

    void set_parameter(const std::string &name, double value) {
        *parameters.addresses()[parameter_index(name)] = value;
    }

    std::map<std::string, double> parameter_map() const {
        std::map<std::string, double> out;
        auto addresses = parameters.addresses();
        for (std::size_t i = 0; i < Params<double>::size; ++i)
            out[Params<double>::names[i]] = *addresses[i];
        return out;
    }

    void set_parameter_map(const std::map<std::string, double> &values) {
        for (const auto &[name, value] : values)
            parameter_index(name);
        for (const auto &[name, value] : values)
            set_parameter(name, value);
    }

protected:
    const Derived &derived() const { return static_cast<const Derived &>(*this); }

#if IMAGINE_HAS_AUTODIFF
    template <int N, typename Evaluate> Eigen::MatrixXd jacobian(Evaluate &&evaluate) const {
        std::vector<std::size_t> columns;
        for (const auto &name : active_parameters)
            columns.push_back(parameter_index(name));
        Params<ad::real> p = parameters.template cast<ad::real>();
        auto addresses = p.addresses();
        Eigen::MatrixXd out(N, columns.size());
        for (std::size_t c = 0; c < columns.size(); ++c) {
            (*addresses[columns[c]])[1] = 1.;
            auto value = evaluate(p);
            for (int k = 0; k < N; ++k)
                out(k, c) = component(value, k)[1];
            (*addresses[columns[c]])[1] = 0.;
        }
        return out;
    }

private:
    static const ad::real &component(const ad::real &value, int) { return value; }
    static const ad::real &component(const std::array<ad::real, 3> &value, int k) { return value[k]; }
#endif
};

template <typename Derived, template <typename> class Params>
class RegularVectorModel : public RegularVectorField, public RegularModel<Derived, Params> {
public:
    Vec3<double> at_position(const double &x, const double &y, const double &z) const override {
        return this->derived().template field<double>(x, y, z, this->parameters);
    }

#if IMAGINE_HAS_AUTODIFF
    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z) const {
        return this->template jacobian<3>(
            [&](const Params<ad::real> &p) { return this->derived().template field<ad::real>(x, y, z, p); });
    }
#endif
};

template <typename Derived, template <typename> class Params>
class RegularScalarModel : public RegularScalarField, public RegularModel<Derived, Params> {
public:
    double at_position(const double &x, const double &y, const double &z) const override {
        return this->derived().template field<double>(x, y, z, this->parameters);
    }

#if IMAGINE_HAS_AUTODIFF
    Eigen::MatrixXd derivative(const double &x, const double &y, const double &z) const {
        return this->template jacobian<1>(
            [&](const Params<ad::real> &p) { return this->derived().template field<ad::real>(x, y, z, p); });
    }
#endif
};

}

#if IMAGINE_HAS_AUTODIFF
#define IMAGINE_INSTANTIATE_AUTODIFF(Model, Result)                     \
    template Result<imagine::ad::real> Model::field<imagine::ad::real>( \
        const double &, const double &, const double &, const Model::parameters_t<imagine::ad::real> &) const;
#else
#define IMAGINE_INSTANTIATE_AUTODIFF(Model, Result)
#endif

#define IMAGINE_INSTANTIATE_MODEL(Model, Result)                                                 \
    template Result<double> Model::field<double>(const double &, const double &, const double &, \
                                                 const Model::parameters_t<double> &) const;     \
    IMAGINE_INSTANTIATE_AUTODIFF(Model, Result)

#define IMAGINE_INSTANTIATE_VECTOR_MODEL(Model) IMAGINE_INSTANTIATE_MODEL(Model, imagine::Vec3)
#define IMAGINE_INSTANTIATE_SCALAR_MODEL(Model) IMAGINE_INSTANTIATE_MODEL(Model, imagine::Scalar)
