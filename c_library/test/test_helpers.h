#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#if autodiff_FOUND
#include <Eigen/Dense>

bool containsNaN(const Eigen::Matrix<autodiff::detail::Real<1, double>, -1, 1, 0, -1, 1>& arr) {
    for (const auto& elem : arr) {
        if (std::isnan(elem.val())) {
            return true;
        }
    }
    return false;
}

#endif


#define assertm(exp, msg) assert(((void)msg, exp))


bool containsNaN(const std::vector<double>& arr) {
    for (const auto& elem : arr) {
        if (std::isnan(elem)) {
            return true;
        }
    }
    return false;
}


bool containsNaN(const std::array<double, 3>& arr) {
    for (const auto& elem : arr) {
        if (std::isnan(elem)) {
            return true;
        }
    }
    return false;
}


// Custom matchers to check array equality

  #if autodiff_FOUND
// Custom matcher for comparing Eigen matrices
template<typename Derived1, typename Derived2>
struct MatrixEqualsMatcher : Catch::Matchers::MatcherBase<Eigen::MatrixBase<Derived1>> {
    
    const Eigen::MatrixBase<Derived2>& m_expected;

    MatrixEqualsMatcher(const Eigen::MatrixBase<Derived2>& expected)
        : m_expected(expected) {}

    // Override the match method to perform the comparison
    bool match(const Eigen::MatrixBase<Derived1>& actual) const override {
        return actual.isApprox(m_expected);
    }

    // Override the describe method to provide a description for the failure message
    std::string describe() const override {
        std::ostringstream ss;
        ss << "is approximately equal to\n" << m_expected;
        return ss.str();
    }

};


// Factory function for creating the custom matcher
template<typename Derived>
MatrixEqualsMatcher<Derived, Derived> EqualsMatrix(const Eigen::MatrixBase<Derived>& expected) {
    return MatrixEqualsMatcher<Derived, Derived>(expected);
}
#endif

template <typename T, size_t N>
// old: struct ArrayEqualsMatcher : Catch::Matchers::Impl::MatcherBase<std::array<T, N>> {
struct ArrayEqualsMatcher : Catch::Matchers::MatcherBase<std::array<T, N>> {
    ArrayEqualsMatcher(const std::array<T, N>& expected) : m_expected(expected) {}

    bool match(const std::array<T, N>& arr) const override {
        for (size_t i = 0; i < N; ++i) {
            if (arr[i] != m_expected[i]) {
                return false;
            }
        }
        return true;
    }

    std::string describe() const override {
        std::ostringstream oss;
        oss << "equals: [";
        for (size_t i = 0; i < N; ++i) {
            if (i > 0) oss << ", ";
            oss << m_expected[i];
        }
        oss << "]";
        return oss.str();
    }

    const std::array<T, N>& m_expected;
};

template <typename T>
// old: struct ArrayEqualsMatcher : Catch::Matchers::Impl::MatcherBase<std::array<T, N>> {
struct PointerEqualsMatcher : Catch::Matchers::MatcherBase<T> {
    PointerEqualsMatcher(const T expected, const size_t N_internal) : m_expected(expected), N2(N_internal) {}

    bool match(const T& arr) const override {
        for (size_t i = 0; i < N2; ++i) {
            if (arr[i] != m_expected[i]) {
                return false;
            }
        }
        return true;
    }

    std::string describe() const override {
        std::ostringstream oss;
        oss << "should equal : [";
        for (size_t i = 0; i < N2; ++i) {
            if (i!=0) oss << ", ";
            oss << m_expected[i];
        }
        oss << "]";
        return oss.str();
    }

    const T m_expected;
    const size_t N2;
};


template <typename T, size_t N>
// old: struct ArrayEqualsMatcher : Catch::Matchers::Impl::MatcherBase<std::array<T, N>> {
struct PointerArrayEqualsMatcher : Catch::Matchers::MatcherBase<std::array<T, N>> {
    PointerArrayEqualsMatcher(const std::array<T, N>& expected, const size_t N_internal) : m_expected(expected), N2(N_internal) {}

    bool match(const std::array<T, N>& arr) const override {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N2; ++j) {
                if (arr[i][j] != m_expected[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    std::string describe() const override {
        std::ostringstream oss;
        oss << "should equal values at : [";
        for (size_t i = 0; i < N; ++i) {
            if (i!=0) oss << ", ";
            oss << m_expected[i];
        }
        oss << "]";
        return oss.str();
    }

    const std::array<T, N>& m_expected;
    const size_t N2;
};

struct VectorEqualsMatcher : Catch::Matchers::MatcherGenericBase {
    VectorEqualsMatcher(const vector& expected) : m_expected(expected), N(m_expected.size()) {}

    bool match(const vector& arr) const {
        for (size_t i = 0; i < N; ++i) {
            if (arr[i] != m_expected[i]) {
                return false;
            }
        }
        return true;
    }

    std::string describe() const override {
        std::ostringstream oss;
        oss << "equals: [";
        for (size_t i = 0; i < N; ++i) {
            if (i > 0) oss << ", ";
            oss << m_expected[i];
        }
        oss << "]";
        return oss.str();
    }

    vector const& m_expected;
    size_t N;
};


// Helper function to create the matcher
template <typename T, size_t N>
ArrayEqualsMatcher<T, N> EqualsArray(const std::array<T, N>& expected) {
    return ArrayEqualsMatcher<T, N>(expected);
}

// Helper function to create the matcher
template <typename T, size_t N>
PointerArrayEqualsMatcher<T, N> EqualsPointerArray(const std::array<T, N>& expected, const size_t N_internal) {
    return PointerArrayEqualsMatcher<T, N>(expected, N_internal);
}

// Helper function to create the matcher
template <typename T>
PointerEqualsMatcher<T> EqualsPointer(const T expected, const size_t N_internal) {
    return PointerEqualsMatcher<T>(expected, N_internal);
}


// Helper function to create the matcher
auto EqualsVector(const vector& expected) -> VectorEqualsMatcher {
    return VectorEqualsMatcher(expected);
}


void _check_array_equality_from_pointer(std::array<double*, 3> a, std::array<double*, 3> b , size_t &n) {
    for (int d = 0; d==3; ++d) {
        std::vector<double>  arr_a(a[d], a[d] + n);
        std::vector<double>  arr_b(b[d], b[d] + n);
        assert (arr_a == arr_b); 
    }
}

