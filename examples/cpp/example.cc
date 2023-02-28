
#include "../../c_library/headers/hamunits.h"
#include "../../c_library/headers/Field.h"
#include "../../c_library/headers/RegularField.h"
#include "../../c_library/headers/RegularJF12.h"
#include "../../c_library/headers/Jaffe.h"
#include "../../c_library/headers/Helix.h"
#include "../../c_library/headers/RandomJF12.h"
#include "../../c_library/headers/EnsslinSteininger.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

void print_pos(std::map <std::string, std::array<double, 3>> pd,
               std::map <std::string, std::shared_ptr<RegularVectorField>> md) {

      auto model_iter = md.begin();


      while (model_iter != md.end()) {
        auto position_iter = pd.begin();
        std::cout << "The model " << model_iter->first << " is evaluated: \n\n";
        while (position_iter != pd.end()) {
          std::array<double, 3> mval = (*(model_iter->second)).getField(position_iter->second);

          std::cout << "Position: ";
          for (size_t l = 0; l < (position_iter->second).size(); l++) {
                            std::cout << (position_iter->second)[l] << " kpc  ";
                          }
          std::cout << "\n";
          std::cout << "Evaluation: ";
          for (size_t l = 0; l < mval.size(); l++) {
                            std::cout << mval[l] << "  ";
                            }
          std::cout << "\n\n";

          ++position_iter;
          }
        std::cout << "\n";
        ++model_iter;
      }
    }

void print_ev_grid_no_grid(std::map <std::string, std::shared_ptr<RegularVectorField>> md, std::array<double, 3> siz) {

  auto model_iter = md.begin();


  while (model_iter != md.end()) {
    std::cout << "The model " << model_iter->first << " is evaluated: \n\n";
    std::array<double*, 3> b_grid = (*(model_iter->second)).on_grid();
    std::cout<< "b_grid: "  <<std::endl;
    int sx = siz[0];
    int sy = siz[1];
    int sz = siz[2];
    for (size_t i = 0; i < sx; i++) {
          for (size_t j = 0; j < sy; j++) {
              for (size_t k = 0; k < sz; k++) {
                  for (size_t l = 0; l < 3; l++) {
                      std::cout << (b_grid[l])[sz*sy*i + sz*j + k] << " ";
                      }
                  std::cout<<"\n";
                  }
              }
          }
      std::cout << "\n";
      ++model_iter;
    }
  }


  void print_ev_grid_irreg(std::map <std::string, std::shared_ptr<RegularVectorField>> md,
  std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) {

  auto model_iter = md.begin();


  while (model_iter != md.end()) {
    std::cout << "The model " << model_iter->first << " is evaluated: \n\n";
    std::array<double*, 3> b_grid = (*(model_iter->second)).on_grid(grid_x, grid_y, grid_z);
    std::cout << "The ref to b_grid in example.cc" << &b_grid << " \n\n";
    std::cout<< "b_grid: "  <<std::endl;
    int sx = grid_x.size();
    int sy = grid_y.size();
    int sz = grid_z.size();
    for (size_t i = 0; i < sx; i++) {
          std::cout<< "b_grid x size: " << sx <<std::endl;
          for (size_t j = 0; j < sy; j++) {
              std::cout<< "b_grid y size: " << sy <<std::endl;
              for (size_t k = 0; k < sz; k++) {
                  std::cout<< "b_grid z size: " << sz <<std::endl;
                  std::cout << "x: "<< grid_x[i] << ", " << "y: "<< grid_y[j]<< ", "<< "z: "<< grid_z[k] << ", ";
                  for (size_t l = 0; l < 3; l++) {
                      std::cout << (b_grid[l])[sz*sy*i + sz*j + k] << " ";
                      }
                  std::cout<<"\n";
                  }
              }
          }
      std::cout << "\n";
      ++model_iter;
      delete b_grid[0];
      delete b_grid[1];
      delete b_grid[2];
    }
  }


  void print_ev_grid_reg(std::map <std::string, std::shared_ptr<RegularVectorField>> md,
  std::array<int, 3> size, std::array<double, 3> zeropoint, std::array<double, 3> increment) {

  auto model_iter = md.begin();


  while (model_iter != md.end()) {
  std::cout << "The model " << model_iter->first << " is evaluated: \n\n";
  std::array<double*, 3> b_grid = (*(model_iter->second)).on_grid(size, zeropoint, increment);
  std::cout<< "b_grid: "  <<std::endl;
  int sx = size[0];
  int sy = size[1];
  int sz = size[2];
  for (size_t i = 0; i < sx; i++) {
        std::cout<< "b_grid x size: " << sx <<std::endl;
        for (size_t j = 0; j < sy; j++) {
            std::cout<< "b_grid y size: " << sy <<std::endl;
            for (size_t k = 0; k < sz; k++) {
                std::cout<< "b_grid z size: " << sz <<std::endl;
                for (size_t l = 0; l < 3; l++) {
                    std::cout << (b_grid[l])[sz*sy*i + sz*j + k] << " ";
                    }
                std::cout<<"\n";
                }
            }
        }
      std::cout << "\n";
      ++model_iter;
      delete b_grid[0];
      delete b_grid[1];
      delete b_grid[2];
    }
  }

int main() {

  std::map <std::string, std::array<double, 3>> position_dict;


  // Define some positions in Galactic cartesian coordinates (units are kpc)
  position_dict["origin"] =  {0., 0., 0.};
  position_dict["pos_on_x_axis"] = {1.7, 0., 0.};
  position_dict["pos_on_y_axis"] = {0., -2.1, 0.};
  position_dict["pos_on_z_axis"] = {0., 0., 1.35};
  position_dict["pos_on_xy_plane"] = {1.9, -2.9, 0.};
  position_dict["pos_on_yz_plane"] = {0., 2.1, -1.6};
  position_dict["pos_on_xz_plane"] = {-2.7, 0., 1.25};

  // Define a irregular grid in Galactic cartesian coordinates (units are kpc)
  const std::vector<double> grid_x {{2., 4., 0., 1., .4}};
  const std::vector<double> grid_y {{4., 6., 0.1, 0., .2}};
  const std::vector<double> grid_z {{-0.2, 0.8, 0.2, 0., 1.}};


  // Define a regular grid in Galactic cartesian coordinates (units are kpc)
  const std::array<int, 3> shape {{4, 3, 2}};
  const std::array<double, 3> increment {{2.1, 0.3, 1.}};
  const std::array<double, 3> zeropoint {{-4., 0.1, -0.3}};

  // Dictionaries
  // using pointer here since RegularField is abstract
  std::map <std::string, std::shared_ptr<RegularVectorField>> reg_mods;

  reg_mods["Jansson Farrar regular"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField());
  reg_mods["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField());
  reg_mods["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField());


  std::map <std::string, std::shared_ptr<RegularVectorField>> reg_mods_w_irreg_grid;

  reg_mods["Jansson Farrar regular"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(grid_x, grid_y, grid_z));
  reg_mods["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(grid_x, grid_y, grid_z));
  reg_mods["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField(grid_x, grid_y, grid_z));


  std::map <std::string, std::shared_ptr<RegularVectorField>> reg_mods_w_reg_grid;

  reg_mods["Jansson Farrar regular"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(shape, zeropoint, increment));
  reg_mods["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(shape, zeropoint, increment));
  reg_mods["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField(shape, zeropoint, increment));


  print_pos(position_dict, reg_mods);


  print_ev_grid_no_grid(reg_mods, {5, 5, 5});

  print_ev_grid_irreg(reg_mods, grid_x, grid_y, grid_z);

  std::map <std::string, std::shared_ptr<RandomField<std::vector<double>, std::vector<double>>>> rand_mods;

  JF12RandomField jf12random();
  ESField esrandom();

  //std::vector<double> jf12_grid = jf12_1.evaluate_model_on_grid(grid_x, grid_y, grid_z);

  //std::vector<std::vector<std::vector<double>>> ymw_grid = ymw16.evaluate_grid(grid_x, grid_y, grid_z);

  //print_ev_grid(jf12_grid, grid_x, grid_y, grid_z);


//  std::cout << "\n b_arm_1: " << jf12.b_arm_1 << std::endl;



  //jf12.b_arm_1 = 40.;

  //std::cout << "\n updated b_arm_1: " << jf12.b_arm_1 << std::endl;

  // Evaluate model again
  //std::vector<double> jf12_val2 = jf12.getField(position_dict.at("pos_on_xy_plane"));
  return 0;
}
