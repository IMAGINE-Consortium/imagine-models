#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "RegularModels.h"

#if autodiff_FOUND

void test_derivative(std::map<std::string, std::map<std::array<double, 3>, Eigen::MatrixXd>> val_pos_map,
                      std::map <std::string, std::shared_ptr<HelixMagneticField>> model_dict
                      ) {
    auto model_iter = model_dict.begin();


    while (model_iter != model_dict.end()) {
        std::map<std::array<double, 3>, Eigen::MatrixXd> val_pos_map_this_model = val_pos_map[(model_iter->first)];
        auto val_pos_iter = val_pos_map_this_model.begin();
        while (val_pos_iter != val_pos_map_this_model.end()) {
          double x = (val_pos_iter->first)[0];
          double y = (val_pos_iter->first)[1];
          double z = (val_pos_iter->first)[2];

          Eigen::MatrixXd mval = (*(model_iter->second)).derivative(x, y, z);
          
          /*std::string assertion_msg1 = "Assert failed with model " + model_iter->first + " and at position " + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z);
          std::string assertion_msg2 = "Assert failed with model eval " + std::to_string(mval[0]) + " " + std::to_string(mval[1]) + " " + std::to_string(mval[2]);
          std::string assertion_msg3 = "But should be " + std::to_string((val_pos_iter->second)[0]) + " " + std::to_string((val_pos_iter->second)[1]) + " " + std::to_string((val_pos_iter->second)[2]);
          std::cout << assertion_msg1  << std::endl;
          std::cout << assertion_msg2  << std::endl;
          std::cout << assertion_msg3  << std::endl;*/
          std::cout << "mval.cols: " << mval.cols()  << " mval.rows: " << mval.rows() << std::endl;
          std::cout << "(val_pos_iter->second)).cols: " <<  (val_pos_iter->second).cols()  << " (val_pos_iter->second).rows: " << (val_pos_iter->second).rows() << std::endl;
          assert(mval == (val_pos_iter->second));
          ++val_pos_iter;
          }
        ++model_iter;
    }
}

int main() {
        std::map<std::string, std::map<std::array<double, 3>, Eigen::MatrixXd>> val_pos_map;

    // Define some positions in Galactic cartesian coordinates (units are kpc)

    std::map<std::array<double, 3>, Eigen::MatrixXd> jf_12_map;
    Eigen::MatrixXd nv{{0., 0., 0.}, 
                        {0., 0., 0.},
                        {0., 0., 0.}};
    jf_12_map[{0., 0., 0.}] = nv; // Galactic center
    jf_12_map[{.1, .3, .4}] = nv; // Within inner boundary
    jf_12_map[{0., 0., 0.}] = nv; // outside outer boundary 


    //std::map<std::array<double, 3>, Eigen::MatrixXd> jaffe_map;
    //jaffe_map[{0., 0., 0.}] = {0., 0., 0.}; // Galactic center

    //std::map<std::array<double, 3>, Eigen::MatrixXd> helix_map;

    val_pos_map["JF12"] = jf_12_map;
    //val_pos_map["Jaffe"] = jaffe_map;
    //val_pos_map["Helix"] = helix_map;


    std::map <std::string, std::shared_ptr<HelixMagneticField>> models;
    models["JF12"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField());
    //models["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField());
    //models["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField());

    test_derivative(val_pos_map, models);

}

#endif