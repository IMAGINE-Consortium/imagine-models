#ifndef REGULARFIELD_H
#define REGULARFIELD_H

#include <vector>
#include <array>
#include <stdexcept>
#include <functional>
#include <iostream>
#include <algorithm>
#include <set>

#include "exceptions.h"
#include "Field.h"

class RegularScalarField : public Field<number, double *>
{
protected:
  // Fields

  // RegularField() : Field<T>() {};
  // Methods
  double *allocate_memory(std::array<int, 3> shp)
  {
    size_t arr_sz = grid_size(shp);
    double *grid_eval = new double[arr_sz];
    return grid_eval;
  }

  void free_memory(double *grid_eval)
  {
    delete grid_eval;
  }

public:
  // -----CONSTRUCTORS-----

  ~RegularScalarField(){};

  RegularScalarField() : Field<number, double *>(){};

  RegularScalarField(std::array<int, 3> shape, std::array<double, 3> reference_point, std::array<double, 3> increment) : Field<number, double *>(shape, reference_point, increment){};

  RegularScalarField(std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : Field<number, double *>(grid_x, grid_y, grid_z){};

  // Fields
  const int ndim = 1;
#if autodiff_FOUND
  const std::set<std::string> all_diff;
  std::set<std::string> active_diff;
#endif
  // Methods

  double *on_grid(int seed = 0)
  {
    if (not initialized_with_grid)
    {
      throw GridException();
    }
    double *grid_eval = allocate_memory(internal_shape);
    if (regular_grid)
    {
      evaluate_function_on_grid<number, double*>(grid_eval, internal_shape, internal_ref_point, internal_increment, [this](double xx, double yy, double zz)
                                        { return at_position(xx, yy, zz); });
    }
    else
    {
      evaluate_function_on_grid<number, double*>(grid_eval, internal_grid_x, internal_grid_y, internal_grid_z, [this](double xx, double yy, double zz)
                                        { return at_position(xx, yy, zz); });
    }
    return grid_eval;
  }

  double *on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, int seed = 0)
  {
    std::array<int, 3> shp = {(int)grid_x.size(), (int)grid_y.size(), (int)grid_z.size()};
    double *grid_eval = allocate_memory(shp);
    evaluate_function_on_grid<number, double*>(grid_eval, grid_x, grid_y, grid_z,
                                      [this](double xx, double yy, double zz)
                                      { return at_position(xx, yy, zz); });
    return grid_eval;
  }

  double *on_grid(const std::array<int, 3> &shape, const std::array<double, 3> &reference_point, const std::array<double, 3> &increment, int seed = 0)
  {
    double *grid_eval = allocate_memory(shape);
    evaluate_function_on_grid<number, double*>(grid_eval, shape, reference_point, increment,
                                      [this](double xx, double yy, double zz)
                                      { return at_position(xx, yy, zz); });
    return grid_eval;
  }

#if autodiff_FOUND

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

class RegularVectorField : public Field<vector, std::array<double *, 3>>
{
protected:
  // Fields

  // RegularField() : Field<T>() {};
  // Methods
  std::array<double *, 3> allocate_memory(std::array<int, 3> shp) override
  {
    std::array<double *, 3> grid_eval;
    size_t arr_sz = grid_size(shp);
    grid_eval[0] = new double[arr_sz];
    grid_eval[1] = new double[arr_sz];
    grid_eval[2] = new double[arr_sz];
    return grid_eval;
  }

  void free_memory(std::array<double *, 3> grid_eval) override
  {
    delete grid_eval[0];
    delete grid_eval[1];
    delete grid_eval[2];
  }

public:
  ~RegularVectorField(){};

  // Constructors
  RegularVectorField() : Field<vector, std::array<double *, 3>>(){};

  RegularVectorField(std::array<int, 3> shape, std::array<double, 3> reference_point, std::array<double, 3> grid_increment) : Field<vector, std::array<double *, 3>>(shape, reference_point, grid_increment){};

  RegularVectorField(std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) : Field<vector, std::array<double *, 3>>(grid_x, grid_y, grid_z){};

  // Fields

  const int ndim = 3;

#if autodiff_FOUND
  const std::set<std::string> all_diff;
  std::set<std::string> active_diff;
#endif
  // Methods

  std::array<double *, 3> on_grid(int seed = 0)
  {
    if (not initialized_with_grid)
    {
      throw GridException();
    }
    std::array<double *, 3> grid_eval = allocate_memory(internal_shape);
    ;
    if (regular_grid)
    {
      evaluate_function_on_grid<vector, std::array<double *, 3>>(grid_eval, internal_shape, internal_ref_point, internal_increment, [this](double xx, double yy, double zz)
                                        { return at_position(xx, yy, zz); });
    }
    else
    {
      evaluate_function_on_grid<vector, std::array<double *, 3>>(grid_eval, internal_grid_x, internal_grid_y, internal_grid_z, [this](double xx, double yy, double zz)
                                        { return at_position(xx, yy, zz); });
    }
    return grid_eval;
  }

  std::array<double *, 3> on_grid(const std::vector<double> &grid_x, const std::vector<double> &grid_y, const std::vector<double> &grid_z, int seed = 0)
  {
    std::array<int, 3> shp = {(int)grid_x.size(), (int)grid_y.size(), (int)grid_z.size()};
    std::array<double *, 3> grid_eval = allocate_memory(shp);
    evaluate_function_on_grid<vector, std::array<double *, 3>>(grid_eval, grid_x, grid_y, grid_z,
                                      [this](double xx, double yy, double zz)
                                      { return at_position(xx, yy, zz); });
    return grid_eval;
  }

  std::array<double *, 3> on_grid(const std::array<int, 3> &shape, const std::array<double, 3> &reference_point, const std::array<double, 3> &increment, int seed = 0)
  {
    std::array<double *, 3> grid_eval = allocate_memory(shape);
    evaluate_function_on_grid<vector, std::array<double *, 3>>(grid_eval, shape, reference_point, increment,
                                      [this](double xx, double yy, double zz)
                                      { return at_position(xx, yy, zz); });
    return grid_eval;
  }

#if autodiff_FOUND

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

#endif
