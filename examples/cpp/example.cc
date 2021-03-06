#include "../../c_library/headers/MagneticField.h"
#include "../../c_library/headers/ThermalElectronField.h"
#include <cassert>
#include <iostream>
#include <vector>

void print_magnetic_pos(std::vector<double> mval, std::vector<double> tp, const char* title) {

      std::cout << "Magnetic field at " << title << " (";
      for (size_t l = 0; l < tp.size(); l++) {
                        std::cout << tp[l] << " kpc, ";
                        }
      std::cout << "):\n";
      for (size_t l = 0; l < mval.size(); l++) {
                        std::cout << mval[l]<< ", ";
                        }
      std::cout << "\n";
      }

void print_thermal_pos(double tval, std::vector<double> tp, const char* title) {

      std::cout << "Theraml electron field at " << title << " (";
      for (size_t l = 0; l < tp.size(); l++) {
                        std::cout << tp[l] << " kpc, ";
                        }
      std::cout << "): " << tval << "\n";
      }

void print_ev_grid(std::vector<std::vector<std::vector<std::vector<double>>>>  b_grid,
  std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) {

  std::cout<< "b_grid: "  <<std::endl;
  for (size_t i = 0; i < b_grid.size(); i++) {
        std::cout<< "b_grid x size: " << b_grid.size() <<std::endl;
        for (size_t j = 0; j < b_grid[i].size(); j++) {
            std::cout<< "b_grid y size: " << b_grid[i].size() <<std::endl;
            for (size_t k = 0; k < b_grid[i][j].size(); k++) {
                std::cout<< "b_grid z size: " << b_grid[i][j].size() <<std::endl;
                std::cout << "x: "<< grid_x[i] << ", " << "y: "<< grid_y[j]<< ", "<< "z: "<< grid_z[k] << ", ";
                for (size_t l = 0; l < b_grid[i][j][k].size(); l++) {
                    std::cout << b_grid[i][j][k][l] << " ";
                    }
                std::cout<<"\n";
                }
            }
        }
    }

int main() {

  JF12MagneticField jf12;
  YMW16ThermalElectronField ymw16;

  // Define some positions in Galactic cartesian coordinates (units are kpc)
  std::vector<double> test_pos{{1., 2., 0.}};
  std::vector<double> origin{{0., 0., 0.}};
  std::vector<double> pos_on_x_axis{{.123, 0., 0.}};
  std::vector<double> pos_on_y_axis{{0., .432, 0.}};
  std::vector<double> pos_on_z_axis{{0., 0., 1.2}};
  std::vector<double> pos_on_xy_plane{{.123, -.213, 0.}};
  std::vector<double> pos_on_yz_plane{{0., -.232, -.3}};
  std::vector<double> pos_on_xz_plane{{-.24, 0., 1.2}};


  std::vector<double> jf12_test_pos = jf12.evaluate_at_pos(test_pos);
  print_magnetic_pos(jf12_test_pos, test_pos, "test_pos");
  std::vector<double> jf12_origin = jf12.evaluate_at_pos(origin);
  print_magnetic_pos(jf12_origin , origin , "origin");
  std::vector<double> jf12_pos_on_x_axis = jf12.evaluate_at_pos(pos_on_x_axis);
  print_magnetic_pos(jf12_pos_on_x_axis, pos_on_x_axis , "pos_on_x_axis");
  std::vector<double> jf12_pos_on_y_axis = jf12.evaluate_at_pos(pos_on_y_axis);
  print_magnetic_pos(jf12_pos_on_y_axis, pos_on_y_axis , "pos_on_y_axis");
  std::vector<double> jf12_pos_on_z_axis = jf12.evaluate_at_pos(pos_on_z_axis);
  print_magnetic_pos(jf12_pos_on_z_axis, pos_on_z_axis , "pos_on_z_axis");
  std::vector<double> jf12_pos_on_xy_plane = jf12.evaluate_at_pos(pos_on_xy_plane);
  print_magnetic_pos(jf12_pos_on_xy_plane, pos_on_xy_plane , "pos_on_xy_plane");
  std::vector<double> jf12_pos_on_yz_plane = jf12.evaluate_at_pos(pos_on_yz_plane);
  print_magnetic_pos(jf12_pos_on_yz_plane, pos_on_yz_plane , "pos_on_yz_plane");
  std::vector<double> jf12_pos_on_xz_plane = jf12.evaluate_at_pos(pos_on_xz_plane);
  print_magnetic_pos(jf12_pos_on_xz_plane, pos_on_xz_plane , "pos_on_xz_plane");


  double ymw_test_pos = ymw16.evaluate_at_pos(test_pos);
  print_thermal_pos(ymw_test_pos, test_pos, "test_pos");
  double ymw_origin = ymw16.evaluate_at_pos(origin);
  print_thermal_pos(ymw_origin , origin , "origin");
  double ymw_pos_on_x_axis = ymw16.evaluate_at_pos(pos_on_x_axis);
  print_thermal_pos(ymw_pos_on_x_axis, pos_on_x_axis , "pos_on_x_axis");
  double ymw_pos_on_y_axis = ymw16.evaluate_at_pos(pos_on_y_axis);
  print_thermal_pos(ymw_pos_on_y_axis, pos_on_y_axis , "pos_on_y_axis");
  double ymw_pos_on_z_axis = ymw16.evaluate_at_pos(pos_on_z_axis);
  print_thermal_pos(ymw_pos_on_z_axis, pos_on_z_axis , "pos_on_z_axis");
  double ymw_pos_on_xy_plane = ymw16.evaluate_at_pos(pos_on_xy_plane);
  print_thermal_pos(ymw_pos_on_xy_plane, pos_on_xy_plane , "pos_on_xy_plane");
  double ymw_pos_on_yz_plane = ymw16.evaluate_at_pos(pos_on_yz_plane);
  print_thermal_pos(ymw_pos_on_yz_plane, pos_on_yz_plane , "pos_on_yz_plane");
  double ymw_pos_on_xz_plane = ymw16.evaluate_at_pos(pos_on_xz_plane);
  print_thermal_pos(ymw_pos_on_xz_plane, pos_on_xz_plane , "pos_on_xz_plane");

  const std::vector<double> grid_x {{2., 4., 0., 1., .4}};
  const std::vector<double> grid_y {{4., 6., 0.1, 0., .2}};
  const std::vector<double> grid_z {{-0.2, 0.8, 0.2, 0., 1.}};


  std::vector<std::vector<std::vector<std::vector<double>>>> jf12_grid = jf12.evaluate_grid(grid_x, grid_y, grid_z);

  std::vector<std::vector<std::vector<double>>> ymw_grid = ymw16.evaluate_grid(grid_x, grid_y, grid_z);

  print_ev_grid(jf12_grid, grid_x, grid_y, grid_z);


  std::cout << "\n b_arm_1: " << jf12.b_arm_1 << std::endl;



  jf12.b_arm_1 = 40.;

  std::cout << "\n updated b_arm_1: " << jf12.b_arm_1 << std::endl;

  // Evaluate moel again
  std::vector<double> jf12_val2 = jf12.evaluate_at_pos(test_pos);
  print_magnetic_pos(jf12_val2, test_pos, "updated jf12, test_pos");
  return 0;
}
