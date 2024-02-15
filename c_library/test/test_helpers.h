#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>


#define assertm(exp, msg) assert(((void)msg, exp))

// Custom matchers to check array equality

struct MatrixEqualsMatcher : Catch::Matchers::MatcherGenericBase {
    MatrixEqualsMatcher(const Eigen::MatrixXd & expected) : m_expected(expected) {}

    bool match(const Eigen::MatrixXd & arr) const {
        
        int r = static_cast<int>(arr.rows());
        if (r != r_expected) return false;
        int c = static_cast<int>(arr.cols());
        if (c != c_expected) return false;
        
        for (size_t i = 0; i < r; ++i) {
            for (size_t j = 0; i < c; ++j) {
                if (arr(i, j) != m_expected(i, j)) {
                    return false;
                }
            }
        }
        
        return true;
    }

    std::string describe() const override {
        std::ostringstream oss;
        oss << "equals: [";
        for (size_t i = 0; i < r_expected; ++i) {
            if (i > 0) oss << "\n ";
            for (size_t j = 0; i < c_expected; ++j) {
                if (j > 0) oss << ", ";
                oss << m_expected(i, j);
            }
        }
        oss << "]";
        return oss.str();
    }

    const Eigen::MatrixXd& m_expected;
    int r_expected = static_cast<int>(m_expected.rows());
    int c_expected = static_cast<int>(m_expected.cols());
};

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
        oss << "equals: [";
        for (size_t i = 0; i < N; ++i) {
            if (i > 0) oss << "\n ";
            for (size_t j = 0; j < N2; ++j) {
                if (j > 0) oss << ", ";
                oss << m_expected[i][j];
            }
        }
        oss << "]";
        return oss.str();
    }

    const std::array<T, N>& m_expected;
    const size_t N2;
};

// old: struct ArrayEqualsMatcher : Catch::Matchers::Impl::MatcherBase<std::array<T, N>> {
struct VectorEqualsMatcher : Catch::Matchers::MatcherGenericBase {
    VectorEqualsMatcher(const vector& expected) : m_expected(expected) {}

    size_t N = m_expected.size();

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
auto EqualsMatrix(const Eigen::MatrixXd& expected) -> MatrixEqualsMatcher {
    return MatrixEqualsMatcher(expected);
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

