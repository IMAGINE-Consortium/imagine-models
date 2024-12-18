#include "RegularModels.h"
//#include <ImagineModelsRandom/RandomModels.h>

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
          std::array<double, 3> mval = (*(model_iter->second)).at_position((position_iter->second)[0], (position_iter->second)[1],(position_iter->second)[2]);

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


template<typename VECTORFIELD>
void print_ev_grid_no_grid(std::map <std::string, std::shared_ptr<VECTORFIELD>> md, std::array<double, 3> siz) {

  auto model_iter = md.begin();


  while (model_iter != md.end()) {
    std::cout << "The model " << model_iter->first << " is evaluated: \n\n";
    std::array<double*, 3> b_grid = (*(model_iter->second)).on_grid(13);
    std::cout<< "b_grid: "  <<std::endl;
    int sx = siz[0];
    int sy = siz[1];
    int sz = siz[2];
    for (size_t i = 0; i < sx; i++) {
          for (size_t j = 0; j < sy; j++) {
              for (size_t k = 0; k < sz; k++) {
                  for (size_t l = 0; l < 3; l++) {
                      std::cout <<  i << ", " <<  j << ", "<<  k << ": " <<  (b_grid[l])[sz*sy*i + sz*j + k] <<  std::endl;
                      }
                  std::cout<<"\n";
                  }
              }
          }
      std::cout << "\n";
      ++model_iter;
    }
    //delete b_grid[0];
    //delete b_grid[1];
    //delete b_grid[2];
  }

template void print_ev_grid_no_grid<RegularVectorField>(std::map <std::string, std::shared_ptr<RegularVectorField>> md, std::array<double, 3> siz); 

//template void print_ev_grid_no_grid<RandomVectorField>(std::map <std::string, std::shared_ptr<RandomVectorField>> md, std::array<double, 3> siz); 

template<typename VECTORFIELD>
void print_ev_grid_reg(std::map <std::string, std::shared_ptr<VECTORFIELD>> md, std::array<int, 3> size, std::array<double, 3> zeropoint, std::array<double, 3> increment) {

  auto model_iter = md.begin();


  while (model_iter != md.end()) {
  std::cout << "The model " << model_iter->first << " is evaluated: \n\n";
  std::array<double*, 3> b_grid = (*(model_iter->second)).on_grid(size, zeropoint, increment, 13);
  std::cout<< "b_grid: "  <<std::endl;
  int sx = size[0];
  int sy = size[1];
  int sz = size[2];
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

template void print_ev_grid_reg<RegularVectorField>(std::map <std::string, std::shared_ptr<RegularVectorField>> md, std::array<int, 3> size, std::array<double, 3> zeropoint, std::array<double, 3> increment);


//template void print_ev_grid_reg<RandomVectorField>(std::map <std::string, std::shared_ptr<RandomVectorField>> md, std::array<int, 3> size, std::array<double, 3> zeropoint, std::array<double, 3> increment);

void print_ev_grid_irreg(std::map <std::string, std::shared_ptr<RegularVectorField>> md, std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z) {

  auto model_iter = md.begin();


  while (model_iter != md.end()) {
    std::cout << "The model " << model_iter->first << " is evaluated: \n\n";
    std::array<double*, 3> b_grid = (*(model_iter->second)).on_grid(grid_x, grid_y, grid_z);
    std::cout<< "b_grid: "  <<std::endl;
    int sx = grid_x.size();
    int sy = grid_y.size();
    int sz = grid_z.size();
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

  reg_mods_w_irreg_grid["Jansson Farrar regular"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(grid_x, grid_y, grid_z));
  reg_mods_w_irreg_grid["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(grid_x, grid_y, grid_z));
  reg_mods_w_irreg_grid["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField(grid_x, grid_y, grid_z));


  std::map <std::string, std::shared_ptr<RegularVectorField>> reg_mods_w_reg_grid;

  reg_mods_w_reg_grid["Jansson Farrar regular"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(shape, zeropoint, increment));
  reg_mods_w_reg_grid["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(shape, zeropoint, increment));
  reg_mods_w_reg_grid["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField(shape, zeropoint, increment));

  //std::map <std::string, std::shared_ptr<RandomVectorField>> rand_mods;

  //rand_mods["Jansson Farrar random"] = std::shared_ptr<JF12RandomField> (new JF12RandomField());
  //rand_mods["Ensslin Steininger"] = std::shared_ptr<ESRandomField> (new ESRandomField());


  //std::map <std::string, std::shared_ptr<RandomVectorField>> rand_mods_w_reg_grid;
  //std::map <std::string, std::shared_ptr<RandomScalarField>> rand_mods_w_reg_grid_scalar;

  //std::cout << "Random scalar fields" << std::endl;
  //rand_mods_w_reg_grid_scalar["Gauss"] = std::shared_ptr<GaussianScalarField> (new GaussianScalarField());

  //std::cout << "Random vector fields" << std::endl;
  //rand_mods_w_reg_grid["Jansson Farrar random"] = std::shared_ptr<JF12RandomField> (new JF12RandomField(shape, zeropoint, increment));
  //rand_mods_w_reg_grid["Ensslin Steininger"] = std::shared_ptr<ESRandomField> (new ESRandomField(shape, zeropoint, increment));



  print_pos(position_dict, reg_mods);


  print_ev_grid_no_grid(reg_mods_w_irreg_grid, {5, 5, 5});




  //rand_mods["Jansson Farrar random"] = std::shared_ptr<JF12RandomField> (new JF12RandomField());
  //rand_mods["Ensslin Steininger"] = std::shared_ptr<ESRandomField> (new ESRandomField());

  //print_ev_grid_reg(rand_mods, shape, zeropoint, increment);
  //print_ev_grid_no_grid(rand_mods_w_reg_grid, {4, 3, 2});

  //print_ev_grid_irreg(reg_mods, grid_x, grid_y, grid_z);



  //std::vector<double> jf12_grid = jf12.on_grid(grid_x, grid_y, grid_z);

  //std::vector<std::vector<std::vector<double>>> ymw_grid = ymw16.evaluate_grid(grid_x, grid_y, grid_z);

  //print_ev_grid_irreg(jf12_grid, grid_x, grid_y, grid_z);


//  std::cout << "\n b_arm_1: " << jf12.b_arm_1 << std::endl;



  //jf12.b_arm_1 = 40.;

  //std::cout << "\n updated b_arm_1: " << jf12.b_arm_1 << std::endl;

  // Evaluate model again
  //std::vector<double> jf12_val2 = jf12.getField(position_dict.at("pos_on_xy_plane"));

  //HelixMagneticField helix = HelixMagneticField();
  //std::cout << helix._derivative(1.9, -2.9, 0., helix.param) << std::endl;
  return 0;
}
