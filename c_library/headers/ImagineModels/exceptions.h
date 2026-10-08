#pragma once

#include <stdexcept>
#include <string>

namespace imagine {

class GridException : public std::invalid_argument {
public:
    GridException(const std::string &message) : std::invalid_argument{message} {}
};

}
