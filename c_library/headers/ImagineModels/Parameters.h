#pragma once

#include <array>
#include <cstddef>

#define IMAGINE_PARAMETER_MEMBER(name, value) T name = value;
#define IMAGINE_PARAMETER_NAME(name, value) #name,
#define IMAGINE_PARAMETER_ADDRESS(name, value) &name,
#define IMAGINE_PARAMETER_CAST(name, value) out.name = U(name);

#define IMAGINE_PARAMETERS(Name, LIST)                                                              \
    template <typename T> struct Name {                                                             \
        LIST(IMAGINE_PARAMETER_MEMBER)                                                              \
        static constexpr std::array names{LIST(IMAGINE_PARAMETER_NAME)};                            \
        static constexpr std::size_t size = names.size();                                           \
        std::array<T *, size> addresses() { return {LIST(IMAGINE_PARAMETER_ADDRESS)}; }             \
        std::array<const T *, size> addresses() const { return {LIST(IMAGINE_PARAMETER_ADDRESS)}; } \
        template <typename U> Name<U> cast() const {                                                \
            Name<U> out;                                                                            \
            LIST(IMAGINE_PARAMETER_CAST)                                                            \
            return out;                                                                             \
        }                                                                                           \
    };
