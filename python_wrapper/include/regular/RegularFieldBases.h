#ifndef REGULARFIELDBASES_H
#define REGULARFIELDBASES_H

#include <pybind11/pybind11.h>

#include "../regular_trampoline.h"
#include "../array_converters.h"

namespace py = pybind11;
using namespace pybind11::literals;

void RegularFieldBases(py::module_ &m) {
        py::class_<RegularVectorField, Field<vector, std::array<double*, 3>>, PyRegularVectorField>(m, "RegularVectorField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("on_grid", [](RegularVectorField &self, py::array_t<double> &grid_x,  py::array_t<double>  &grid_y, py::array_t<double>  &grid_z)  {
            size_t sx = grid_x.size();
            size_t sy = grid_y.size();
            size_t sz = grid_z.size();
            std::vector<double> grid_x_vec{grid_x.data(), grid_x.data() + sx}; 
            std::vector<double> grid_y_vec{grid_y.data(), grid_y.data() + sy}; 
            std::vector<double> grid_z_vec{grid_z.data(), grid_z.data() + sz}; 
            std::array<double*, 3> f = self.on_grid(grid_x_vec, grid_y_vec, grid_z_vec);
            auto li = from_pointer_array_to_list_pyarray(f, sx, sy, sz);
            //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
            return li;},
            py::arg("grid_x").noconvert(), py::arg("grid_y").noconvert(), py::arg("grid_z").noconvert(), py::return_value_policy::take_ownership)


        .def("on_grid", [](RegularVectorField &self, std::array<int, 3> &shape,  std::array<double, 3>  &reference_point, std::array<double, 3>  &increment)  {
          std::array<double*, 3> f = self.on_grid(shape, reference_point, increment);
          size_t sx = shape[0];
          size_t sy = shape[1];
          size_t sz = shape[2];
          auto arr = from_pointer_array_to_list_pyarray(f, sx, sy, sz);
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;}, 
          py::kw_only(), py::arg("shape").noconvert(), py::arg("reference_point"), py::arg("increment"), py::return_value_policy::take_ownership)

        .def("on_grid", [](RegularVectorField &self)  {
          std::array<double*, 3> f = self.on_grid();
          size_t sx = self.internal_shape[0];
          size_t sy = self.internal_shape[1];
          size_t sz = self.internal_shape[2];
          auto arr = from_pointer_array_to_list_pyarray(f, sx, sy, sz);
          //py::array_t<double> arr = py::array(f.size(), f.data());  // produces a copy!
          return arr;}, py::return_value_policy::take_ownership);
        

// Regular Scalar Base Class
    py::class_<RegularScalarField, Field<number, double*>, PyRegularScalarField>(m, "RegularScalarField")
        .def(py::init<>())
        .def(py::init<std::vector<double> &, std::vector<double> &, std::vector<double> &>())
        .def(py::init<std::array<int, 3> &, std::array<double, 3> &, std::array<double, 3> &>())

        .def("on_grid", [](RegularScalarField &self, py::array_t<double> &grid_x,  py::array_t<double>  &grid_y, py::array_t<double>  &grid_z)  {
          size_t sx = grid_x.size();
          size_t sy = grid_y.size();
          size_t sz = grid_z.size();
          std::vector<double> grid_x_vec{grid_x.data(), grid_x.data() + sx}; 
          std::vector<double> grid_y_vec{grid_y.data(), grid_y.data() + sy}; 
          std::vector<double> grid_z_vec{grid_z.data(), grid_z.data() + sz}; 
          double* f = self.on_grid(grid_x_vec, grid_y_vec, grid_z_vec);
          auto arr = from_pointer_to_pyarray(std::move(f), sx, sy, sz);
          //arr.resize({sx, sy, sz});
          return arr;},
          py::arg("grid_x").noconvert(), py::arg("grid_y").noconvert(), py::arg("grid_z").noconvert(), py::return_value_policy::take_ownership)
        
        .def("on_grid", [](RegularScalarField &self, std::array<int, 3>  shape, std::array<double, 3>  reference_point, std::array<double, 3>  increment)  {
          double* f = self.on_grid(shape, reference_point, increment);
          size_t sx = shape[0];
          size_t sy = shape[1];
          size_t sz = shape[2];
          auto arr = from_pointer_to_pyarray(std::move(f), sx, sy, sz);
          return arr;},
          py::kw_only(), py::arg("shape").noconvert(), py::arg("reference_point"), py::arg("increment"), 
          py::return_value_policy::take_ownership)


        .def("on_grid", [](RegularScalarField &self)  {
          double* f = self.on_grid();
          size_t sx = self.internal_shape[0];
          size_t sy = self.internal_shape[1];
          size_t sz = self.internal_shape[2];
          auto arr = from_pointer_to_pyarray(std::move(f), sx, sy, sz);
          return arr;}, py::return_value_policy::take_ownership);
}
#endif