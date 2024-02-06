#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>

#include "RegularModels.h"

#define assertm(exp, msg) assert(((void)msg, exp))


void _check_array_equality_from_pointer(std::array<double*, 3> a, std::array<double*, 3> b , size_t &n) {
    for (int d = 0; d==3; ++d) {
        std::vector<double>  arr_a(a[d], a[d] + n);
        std::vector<double>  arr_b(b[d], b[d] + n);
        assert (arr_a == arr_b); 
    }

}




void test_grid(std::map <std::string, std::shared_ptr<RegularVectorField>> models_no_grid, 
               std::map <std::string, std::shared_ptr<RegularVectorField>> models_regular_grid, 
               std::map <std::string, std::shared_ptr<RegularVectorField>> models_irregular_grid, 
               std::array<int, 3> shape, std::array<double, 3> refpoint, std::array<double, 3> increment, 
               std::vector<double> grid_x, std::vector<double> grid_y, std::vector<double> grid_z
              ) {
    auto model_iter = models_no_grid.begin();

    


    while (model_iter != models_no_grid.end()) { 
        std::array<double*, 3> eval_irregular = (*models_irregular_grid[model_iter->first]).on_grid(); 
        std::array<double*, 3> eval_regular = (*models_regular_grid[model_iter->first]).on_grid(); 
        std::array<double*, 3> eval_no_grid_regular = (*models_no_grid[model_iter->first]).on_grid(shape, refpoint, increment); 
        std::array<double*, 3> eval_no_grid_irregular = (*models_no_grid[model_iter->first]).on_grid(grid_x, grid_y, grid_z); 
        size_t irreg_n = grid_x.size()*grid_y.size()*grid_z.size();
        size_t reg_n = shape[0]*shape[1]*shape[2];
        // std::array<double*, 3> eval_no_grid = (*models_no_grid[model_iter->first]).on_grid(); should raise an exception

        _check_array_equality_from_pointer(eval_irregular, eval_no_grid_irregular, irreg_n);
        _check_array_equality_from_pointer(eval_regular, eval_no_grid_regular, reg_n);
        ++model_iter;
    }

}

     


int main() {


    // Define a irregular grid in Galactic cartesian coordinates (units are kpc)
    const std::vector<double> grid_x {{2., 4., 0., 1., .4}};
    const std::vector<double> grid_y {{4., 6., 0.1, 0., .2}};
    const std::vector<double> grid_z {{-0.2, 0.8, 0.2, 0., 1.}};


    // Define a regular grid in Galactic cartesian coordinates (units are kpc)
    const std::array<int, 3> shape {{4, 3, 2}};
    const std::array<double, 3> refpoint {{-4., 0.1, -0.3}};
    const std::array<double, 3> increment {{2.1, 0.3, 1.}};


    // Dictionaries
    // using pointer here since RegularField is abstract
    std::map <std::string, std::shared_ptr<RegularVectorField>> models_w_empty_constructor;
    std::map <std::string, std::shared_ptr<RegularVectorField>> models_w_regular_constructor;
    std::map <std::string, std::shared_ptr<RegularVectorField>> models_w_irregular_constructor;

    models_w_empty_constructor["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField());
    //models_w_empty_constructor["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField());
    models_w_empty_constructor["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField());

    models_w_regular_constructor["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(shape, refpoint, increment));
    //models_w_regular_constructor["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField(shape, refpoint, increment));

    models_w_irregular_constructor["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(grid_x, grid_y, grid_z));
    //models_w_irregular_constructor["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["Helix"] = std::shared_ptr<HelixMagneticField> (new HelixMagneticField(grid_x, grid_y, grid_z));



    test_grid(models_w_empty_constructor, models_w_regular_constructor, models_w_irregular_constructor, shape, refpoint, increment, grid_x, grid_y, grid_z);
}



