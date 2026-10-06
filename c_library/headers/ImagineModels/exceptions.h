#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <stdexcept>
#include <string>

namespace imagine {

class GridException : public std::invalid_argument
{
public:
    GridException (const std::string &message) : std::invalid_argument{message} {}
};

class NotImplementedException : public std::logic_error
{
public:
    NotImplementedException () : std::logic_error{"Function not yet implemented."} {}
};


class FieldException : public std::logic_error
{
public:
    FieldException () : std::logic_error{"Field cannot be intialized this way."} {}
};


class DivergenceException : public std::logic_error
{
public:
    DivergenceException () : std::logic_error{"The divergence of a vectorfield can only be calculated in 3 dimensions"} {}
};

}

#endif