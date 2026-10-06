#ifndef REGULARFIELD_H
#define REGULARFIELD_H

#include <cstddef>

#include "ImagineModels/types.h"
#include "ImagineModels/Grid.h"

namespace imagine {

class RegularScalarField
{
public:
  virtual ~RegularScalarField() = default;

  virtual double at_position(const double &x, const double &y, const double &z) const = 0;

  ScalarGridData evaluate(const Grid &grid) const
  {
    ScalarGridData out(grid_shape(grid));
    for_each_point(grid, [&](std::size_t idx, double x, double y, double z)
                   { out(0, idx) = at_position(x, y, z); });
    return out;
  }
};

class RegularVectorField
{
public:
  virtual ~RegularVectorField() = default;

  virtual Vec3<double> at_position(const double &x, const double &y, const double &z) const = 0;

  VectorGridData evaluate(const Grid &grid) const
  {
    VectorGridData out(grid_shape(grid));
    for_each_point(grid, [&](std::size_t idx, double x, double y, double z)
                   {
                     Vec3<double> v = at_position(x, y, z);
                     out(0, idx) = v[0];
                     out(1, idx) = v[1];
                     out(2, idx) = v[2]; });
    return out;
  }
};

}

#endif
