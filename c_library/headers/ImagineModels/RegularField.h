#ifndef REGULARFIELD_H
#define REGULARFIELD_H

#include <vector>
#include <array>
#include <string>
#include <set>

#include "ImagineModels/types.h"
#include "ImagineModels/Grid.h"

namespace imagine {

class RegularScalarField
{
public:
  virtual ~RegularScalarField() = default;

#if IMAGINE_HAS_AUTODIFF
  const std::set<std::string> all_diff;
  std::set<std::string> active_diff;
#endif

  virtual number at_position(const double &x, const double &y, const double &z) const = 0;

  ScalarGridData evaluate(const Grid &grid) const
  {
    ScalarGridData out(grid_shape(grid));
    for_each_point(grid, [&](std::size_t idx, double x, double y, double z)
                   { out(0, idx) = static_cast<double>(at_position(x, y, z)); });
    return out;
  }

#if IMAGINE_HAS_AUTODIFF

  Eigen::VectorXd _filter_diff(Eigen::VectorXd inp) const
  {
    if (active_diff.size() != all_diff.size())
    {
      std::vector<int> i_to_keep;
      for (std::string s : active_diff)
      {
        if (auto search = all_diff.find(s); search != all_diff.end())
        {
          int index = std::distance(all_diff.begin(), search);
          i_to_keep.push_back(index);
        }
      }
      return inp(Eigen::all, i_to_keep);
    }
    return inp;
  }

#endif
};

class RegularVectorField
{
public:
  virtual ~RegularVectorField() = default;

#if IMAGINE_HAS_AUTODIFF
  const std::set<std::string> all_diff;
  std::set<std::string> active_diff;
#endif

  virtual vector at_position(const double &x, const double &y, const double &z) const = 0;

  VectorGridData evaluate(const Grid &grid) const
  {
    VectorGridData out(grid_shape(grid));
    for_each_point(grid, [&](std::size_t idx, double x, double y, double z)
                   {
                     vector v = at_position(x, y, z);
                     out(0, idx) = static_cast<double>(v[0]);
                     out(1, idx) = static_cast<double>(v[1]);
                     out(2, idx) = static_cast<double>(v[2]); });
    return out;
  }

#if IMAGINE_HAS_AUTODIFF

  Eigen::MatrixXd _filter_diff(Eigen::MatrixXd inp) const
  {
    if (active_diff.size() != all_diff.size())
    {
      std::vector<int> i_to_keep;
      for (std::string s : active_diff)
      {
        if (auto search = all_diff.find(s); search != all_diff.end())
        {
          int index = std::distance(all_diff.begin(), search);
          i_to_keep.push_back(index);
        }
      }
      return inp(Eigen::all, i_to_keep);
    }
    return inp;
  }

#endif
};

}

#endif
