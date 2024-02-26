#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>
#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <memory>


#include "RandomModels.h"
#include "../test_helpers.h"


TEST_CASE("RandomVectorFields") { 
    //                                   

    // Define a regular grid in Galactic cartesian coordinates (units are kpc)
    const std::array<int, 3> shape {{4, 3, 2}};
    const std::array<double, 3> refpoint {{-4., 0.1, -0.3}};
    const std::array<double, 3> increment {{2.1, 0.3, 1.}};


    // Dictionaries of all models, uniform an helix are tested elsewhere
    // using pointer here since RegularField is abstract
    std::map <std::string, std::shared_ptr<RandomVectorField>> models_w_empty_constructor;
    std::map <std::string, std::shared_ptr<RandomVectorField>> models_w_regular_constructor;

    models_w_empty_constructor["JF12"] = std::shared_ptr<JF12RandomField> (new JF12RandomField());
    models_w_regular_constructor["JF12"] = std::shared_ptr<JF12RandomField> (new JF12RandomField(shape, refpoint, increment));
	models_w_irregular_constructor["JF12"] = std::shared_ptr<JF12RandomField> (new JF12RandomField(positions_x, positions_y, positions_z));

    auto empty_model_iter = models_w_empty_constructor.begin();
    auto regular_model_iter = models_w_empty_constructor.begin();

    while (empty_model_iter != models_w_empty_constructor.end()) { 

        UNSCOPED_INFO("RandomModels testing: on grid failed for: ");
        CAPTURE(irregular_model_iter->first);
        //std::cout << irregular_model_iter->first << std::endl;
        
        REQUIRE_THROWS_WITH((*models_w_empty_constructor[empty_model_iter->first]).on_grid(), "The class has not been initialized with a grid, hence on_grid can only be called with a grid provided.");
        std::array<double*, 3> eval_irregular = (*models_w_irregular_constructor[irregular_model_iter->first]).
        on_grid(); 
        std::array<double*, 3> eval_regular = (*models_w_regular_constructor[regular_model_iter->first]).on_grid(); 
        std::array<double*, 3> eval_no_grid_regular = (*models_w_empty_constructor[empty_model_iter->first]).on_grid(shape, refpoint, increment); 


        REQUIRE_THAT(eval_regular, EqualsPointerArray(eval_no_grid_regular, 4*3*2));
        
        ++empty_model_iter;
        ++regular_model_iter;
    }
}
