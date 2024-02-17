#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <memory>


#include "RegularModels.h"
#include "../test_helpers.h"


TEST_CASE("RegularMagneticField") { 

    const std::vector<double> grid_x {{2., 4., 0., 1., .4}};
    const std::vector<double> grid_y {{4., 6., 0.1, 0., .2}};
    const std::vector<double> grid_z {{-0.2, 0.8, 0.2, 0., 1.}};


    // Define a regular grid in Galactic cartesian coordinates (units are kpc)
    const std::array<int, 3> shape {{4, 3, 2}};
    const std::array<double, 3> refpoint {{-4., 0.1, -0.3}};
    const std::array<double, 3> increment {{2.1, 0.3, 1.}};


    // Dictionaries of all models, uniform an helix are tested elsewhere
    // using pointer here since RegularField is abstract
    std::map <std::string, std::shared_ptr<RegularVectorField>> models_w_empty_constructor;
    std::map <std::string, std::shared_ptr<RegularVectorField>> models_w_regular_constructor;
    std::map <std::string, std::shared_ptr<RegularVectorField>> models_w_irregular_constructor;

    models_w_empty_constructor["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField());
    models_w_empty_constructor["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField());
    models_w_empty_constructor["Sun2007"] = std::shared_ptr<SunMagneticField> (new SunMagneticField());
    models_w_empty_constructor["Han2018"] = std::shared_ptr<HanMagneticField> (new HanMagneticField());
    models_w_empty_constructor["Pshirkov"] = std::shared_ptr<PshirkovMagneticField> (new PshirkovMagneticField());
    models_w_empty_constructor["HMR"] = std::shared_ptr<HMRMagneticField> (new HMRMagneticField());
    models_w_empty_constructor["TT"] = std::shared_ptr<TTMagneticField> (new TTMagneticField());
    models_w_empty_constructor["TF"] = std::shared_ptr<TFMagneticField> (new TFMagneticField());
    models_w_empty_constructor["Fauvet"] = std::shared_ptr<FauvetMagneticField> (new FauvetMagneticField());
    models_w_empty_constructor["Stanev"] = std::shared_ptr<StanevBSSMagneticField> (new StanevBSSMagneticField());
    models_w_empty_constructor["WMAP"] = std::shared_ptr<WMAPMagneticField> (new WMAPMagneticField());
    models_w_empty_constructor["Archimedean spiral"] = std::shared_ptr<ArchimedeanMagneticField> (new ArchimedeanMagneticField());

    models_w_regular_constructor["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(shape, refpoint, increment));
    models_w_regular_constructor["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["Sun2007"] = std::shared_ptr<SunMagneticField> (new SunMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["Han2018"] = std::shared_ptr<HanMagneticField> (new HanMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["Pshirkov"] = std::shared_ptr<PshirkovMagneticField> (new PshirkovMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["HMR"] = std::shared_ptr<HMRMagneticField> (new HMRMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["TT"] = std::shared_ptr<TTMagneticField> (new TTMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["TF"] = std::shared_ptr<TFMagneticField> (new TFMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["Fauvet"] = std::shared_ptr<FauvetMagneticField> (new FauvetMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["Stanev"] = std::shared_ptr<StanevBSSMagneticField> (new StanevBSSMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["WMAP"] = std::shared_ptr<WMAPMagneticField> (new WMAPMagneticField(shape, refpoint, increment));
    models_w_regular_constructor["Archimedean spiral"] = std::shared_ptr<ArchimedeanMagneticField> (new ArchimedeanMagneticField(shape, refpoint, increment));

	models_w_irregular_constructor["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["Sun2007"] = std::shared_ptr<SunMagneticField> (new SunMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["Han2018"] = std::shared_ptr<HanMagneticField> (new HanMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["Pshirkov"] = std::shared_ptr<PshirkovMagneticField> (new PshirkovMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["HMR"] = std::shared_ptr<HMRMagneticField> (new HMRMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["TT"] = std::shared_ptr<TTMagneticField> (new TTMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["TF"] = std::shared_ptr<TFMagneticField> (new TFMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["Fauvet"] = std::shared_ptr<FauvetMagneticField> (new FauvetMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["Stanev"] = std::shared_ptr<StanevBSSMagneticField> (new StanevBSSMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["WMAP"] = std::shared_ptr<WMAPMagneticField> (new WMAPMagneticField(grid_x, grid_y, grid_z));
    models_w_irregular_constructor["Archimedean spiral"] = std::shared_ptr<ArchimedeanMagneticField> (new ArchimedeanMagneticField(grid_x, grid_y, grid_z));

    auto empty_model_iter = models_w_empty_constructor.begin();
    auto regular_model_iter = models_w_empty_constructor.begin();
    auto irregular_model_iter = models_w_empty_constructor.begin();

    while (empty_model_iter != models_w_empty_constructor.end()) { 
        REQUIRE_THROWS_WITH((*models_w_empty_constructor[empty_model_iter->first]).on_grid(), "The class has not been initialized with a grid, hence on_grid can only be called with a grid provided.");
        
        std::array<double*, 3> eval_irregular = (*models_w_irregular_constructor[irregular_model_iter->first]).
        on_grid(); 
        std::array<double*, 3> eval_regular = (*models_w_regular_constructor[regular_model_iter->first]).on_grid(); 
        std::array<double*, 3> eval_no_grid_regular = (*models_w_empty_constructor[empty_model_iter->first]).on_grid(shape, refpoint, increment); 
        std::array<double*, 3> eval_no_grid_irregular = (*models_w_empty_constructor[empty_model_iter->first]).on_grid(grid_x, grid_y, grid_z); 

        REQUIRE_THAT(eval_regular, EqualsPointerArray(eval_no_grid_regular, 5));
        REQUIRE_THAT(eval_irregular, EqualsPointerArray(eval_no_grid_irregular, 5));
        
        REQUIRE_NOTHROW((*models_w_empty_constructor[empty_model_iter->first]).at_position(1., 2., -3.1));
        REQUIRE_NOTHROW((*models_w_irregular_constructor[irregular_model_iter->first]).at_position(0., 0., 0.));
        REQUIRE_NOTHROW((*models_w_regular_constructor[regular_model_iter->first]).at_position(-10., 23., -33.1));

        #if autodiff_FOUND
        // test derivative
        REQUIRE_NOTHROW((*models_w_empty_constructor[empty_model_iter->first]).derivative(1., 2., -3.1));
        REQUIRE_NOTHROW((*models_w_irregular_constructor[irregular_model_iter->first]).derivative(0., 0., 0.));
        REQUIRE_NOTHROW((*models_w_regular_constructor[regular_model_iter->first]).derivative(-10., 23., -33.1));
        #endif

        ++empty_model_iter;
        ++regular_model_iter;
        ++irregular_model_iter;
    }

}