#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "RegularModels.h"

void test_at_position(std::map<std::string, std::map<std::array<double, 3>, vector>> val_pos_map,
                      std::map <std::string, std::shared_ptr<RegularVectorField>> model_dict
                      ) {
    auto model_iter = model_dict.begin();


    while (model_iter != model_dict.end()) {
        std::map<std::array<double, 3>, vector> val_pos_map_this_model = val_pos_map[(model_iter->first)];
        auto val_pos_iter = val_pos_map_this_model.begin();
        while (val_pos_iter != val_pos_map_this_model.end()) {
          double x = (val_pos_iter->first)[0];
          double y = (val_pos_iter->first)[1];
          double z = (val_pos_iter->first)[2];

          vector mval = (*(model_iter->second)).at_position(x, y, z);
          /*std::string assertion_msg1 = "Assert failed with model " + model_iter->first + " and at position " + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z);
          std::string assertion_msg2 = "Assert failed with model eval " + std::to_string(mval[0]) + " " + std::to_string(mval[1]) + " " + std::to_string(mval[2]);
          std::string assertion_msg3 = "But should be " + std::to_string((val_pos_iter->second)[0]) + " " + std::to_string((val_pos_iter->second)[1]) + " " + std::to_string((val_pos_iter->second)[2]);
          std::cout << assertion_msg1  << std::endl;
          std::cout << assertion_msg2  << std::endl;
          std::cout << assertion_msg3  << std::endl;*/
          assert(mval == (val_pos_iter->second));
          ++val_pos_iter;
          }
        ++model_iter;
    }
}

int main() {
    // Define some positions in Galactic cartesian coordinates (units are kpc)
    vector zv{{0., 0., 0.}};

    vector z1{{0., 0., 1.}};

    std::map<std::array<double, 3>, vector> jf_12_map;
    jf_12_map[{0., 0., 0.}] = zv; // Galactic center
    jf_12_map[{.1, .3, .4}] = zv; // Within inner boundary
    jf_12_map[{0., 0., 0.}] = zv; // outside outer boundary 


    std::map<std::array<double, 3>, vector> jaffe_map;
    jaffe_map[{0., 0., 0.}] = zv; // Galactic center
    jaffe_map[{0., 0., 1.}] = zv; // Galactic center

    std::map<std::array<double, 3>, vector> pshirkov_map;
    pshirkov_map[{0., 0., 0.}] = zv; // Galactic center

    std::map<std::array<double, 3>, vector> helix_map;

    std::map<std::string, std::map<std::array<double, 3>, vector>> val_pos_map;
    val_pos_map["JF12"] = jf_12_map;
    val_pos_map["Jaffe"] = jaffe_map;
    val_pos_map["Helix"] = helix_map;
    val_pos_map["Pshirkov"] = pshirkov_map;


    std::map <std::string, std::shared_ptr<RegularVectorField>> models;
    models["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField());
    models["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField());
    models["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField());
    models["Pshirkov"] = std::shared_ptr<PshirkovMagneticField> (new PshirkovMagneticField());

    test_at_position(val_pos_map, models);

}