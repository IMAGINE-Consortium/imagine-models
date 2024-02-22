#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>
#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <memory>


#include "RegularModels.h"
#include "../test_helpers.h"


TEST_CASE("RegularVectorFields") { 
    //                                   
    const std::vector<double> positions_x {{0., 0. , 0.   , 0.  , 0.  , 0.3, -17.2, 
    0.     , 0.    , 0.   , 0.    , -35.032, 2.2, -6. , .234,  -0.00424, -2.2, -2.32 , .234,  
    0.43 , 8.23    , 3.   , 1., -0.001   , -5.342, -5.23, -2.432  }};
    const std::vector<double> positions_y {{0., 0. , 0.   , 11.5, -7.2, 0. , 0.   , 
    -0.0032, 23.353, -6.1 , 4.    , 0.     , 0. , 0.  , 0.  ,  -0.00132, 2.2 ,  -6.42, .234,   
    6.123, 8.1     , -3.  , -1., 0.0033  , 0.45  , -30.2, -1.543  }};
    const std::vector<double> positions_z {{0., 3.1, -0.01, 0.  , 0.  , 0. , 0.   , 
    -5.32  , 0.012 , 34.02, -8.431, -5.132, 9.3 , -6.1, 4.  , 0.       , 0.  , 0.    , 0.  ,
    10.3 , -0.02   , 3.   , -1., 0.000432, -15.43, 0.001, -7.123  }};

    int size_pos = positions_x.size();

    UNSCOPED_INFO("Test case broken!");
    REQUIRE(size_pos == positions_y.size());
    REQUIRE(size_pos == positions_z.size());

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

	models_w_irregular_constructor["JF12"] = std::shared_ptr<JF12MagneticField> (new JF12MagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["Jaffe"] = std::shared_ptr<JaffeMagneticField> (new JaffeMagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["Sun2007"] = std::shared_ptr<SunMagneticField> (new SunMagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["Han2018"] = std::shared_ptr<HanMagneticField> (new HanMagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["Pshirkov"] = std::shared_ptr<PshirkovMagneticField> (new PshirkovMagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["HMR"] = std::shared_ptr<HMRMagneticField> (new HMRMagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["TT"] = std::shared_ptr<TTMagneticField> (new TTMagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["TF"] = std::shared_ptr<TFMagneticField> (new TFMagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["Fauvet"] = std::shared_ptr<FauvetMagneticField> (new FauvetMagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["Stanev"] = std::shared_ptr<StanevBSSMagneticField> (new StanevBSSMagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["WMAP"] = std::shared_ptr<WMAPMagneticField> (new WMAPMagneticField(positions_x, positions_y, positions_z));
    models_w_irregular_constructor["Archimedean spiral"] = std::shared_ptr<ArchimedeanMagneticField> (new ArchimedeanMagneticField(positions_x, positions_y, positions_z));

    auto empty_model_iter = models_w_empty_constructor.begin();
    auto regular_model_iter = models_w_empty_constructor.begin();
    auto irregular_model_iter = models_w_empty_constructor.begin();

    while (empty_model_iter != models_w_empty_constructor.end()) { 
        auto pos_x_iter = positions_x.begin();
        auto pos_y_iter = positions_y.begin();
        auto pos_z_iter = positions_z.begin();
        while (pos_x_iter != positions_x.end()) { 
            UNSCOPED_INFO("RegularModels testing: at_position call failed for: ");
            CAPTURE(irregular_model_iter->first, *pos_x_iter, *pos_y_iter, *pos_z_iter);
            REQUIRE_NOTHROW((*models_w_empty_constructor[empty_model_iter->first]).at_position(*pos_x_iter, *pos_y_iter, *pos_z_iter));
            UNSCOPED_INFO("RegularModels testing: at_position contained NaN: ");
            CAPTURE(irregular_model_iter->first, *pos_x_iter, *pos_y_iter, *pos_z_iter);
            auto eval = (*models_w_empty_constructor[empty_model_iter->first]).at_position(*pos_x_iter, *pos_y_iter, *pos_z_iter);
            REQUIRE_FALSE(containsNaN(eval));
            
            #if autodiff_FOUND
            // test derivative
            UNSCOPED_INFO("RegularModels testing: derivative call failed for: ");
            CAPTURE(irregular_model_iter->first, pos_x_iter, pos_y_iter, pos_z_iter);
            REQUIRE_NOTHROW((*models_w_empty_constructor[empty_model_iter->first]).derivative(*pos_x_iter, *pos_y_iter, *pos_z_iter));
            #endif

            ++pos_x_iter;
            ++pos_y_iter;
            ++pos_z_iter;
            }

        UNSCOPED_INFO("RegularModels testing: on grid failed for: ");
        CAPTURE(irregular_model_iter->first);
        //std::cout << irregular_model_iter->first << std::endl;
        
        REQUIRE_THROWS_WITH((*models_w_empty_constructor[empty_model_iter->first]).on_grid(), "The class has not been initialized with a grid, hence on_grid can only be called with a grid provided.");
        std::array<double*, 3> eval_irregular = (*models_w_irregular_constructor[irregular_model_iter->first]).
        on_grid(); 
        std::array<double*, 3> eval_regular = (*models_w_regular_constructor[regular_model_iter->first]).on_grid(); 
        std::array<double*, 3> eval_no_grid_regular = (*models_w_empty_constructor[empty_model_iter->first]).on_grid(shape, refpoint, increment); 
        std::array<double*, 3> eval_no_grid_irregular = (*models_w_empty_constructor[empty_model_iter->first]).on_grid(positions_x, positions_y, positions_z);         

        CAPTURE(eval_irregular[0][0], eval_irregular[1][0], eval_irregular[2][0], 
        eval_no_grid_irregular[0][0], eval_no_grid_irregular[1][0], eval_no_grid_irregular[2][0],
        eval_irregular[0][1], eval_irregular[1][1], eval_irregular[2][1], 
        eval_no_grid_irregular[0][1], eval_no_grid_irregular[1][1], eval_no_grid_irregular[2][1],
        eval_irregular[0][2], eval_irregular[1][2], eval_irregular[2][2], 
        eval_no_grid_irregular[0][2], eval_no_grid_irregular[1][2], eval_no_grid_irregular[2][2],
        eval_irregular[0][3], eval_irregular[1][3], eval_irregular[2][3], 
        eval_no_grid_irregular[0][3], eval_no_grid_irregular[1][3], eval_no_grid_irregular[2][3],
        eval_irregular[0][4], eval_irregular[1][4], eval_irregular[2][4], 
        eval_no_grid_irregular[0][4], eval_no_grid_irregular[1][4], eval_no_grid_irregular[2][4]
        );

        REQUIRE_THAT(eval_regular, EqualsPointerArray(eval_no_grid_regular, 4*3*2));
        REQUIRE_THAT(eval_irregular, EqualsPointerArray(eval_no_grid_irregular, size_pos));
        
        ++empty_model_iter;
        ++regular_model_iter;
        ++irregular_model_iter;
    }
}

TEST_CASE("RegularScalarFields") { 
    //                                   
    const std::vector<double> positions_x {{0., 0. , 0.   , 0.  , 0.  , 0.3, -17.2, 
    0.     , 0.    , 0.   , 0.    , -35.032, 2.2, -6. , .234,  -0.00424, -2.2, -2.32 , .234,  
    0.43 , 8.23    , 3.   , 1., -0.001   , -5.342, -5.23, -2.432  }};
    const std::vector<double> positions_y {{0., 0. , 0.   , 11.5, -7.2, 0. , 0.   , 
    -0.0032, 23.353, -6.1 , 4.    , 0.     , 0. , 0.  , 0.  ,  -0.00132, 2.2 ,  -6.42, .234,   
    6.123, 8.1     , -3.  , -1., 0.0033  , 0.45  , -30.2, -1.543  }};
    const std::vector<double> positions_z {{0., 3.1, -0.01, 0.  , 0.  , 0. , 0.   , 
    -5.32  , 0.012 , 34.02, -8.431, -5.132, 9.3 , -6.1, 4.  , 0.       , 0.  , 0.    , 0.  ,
    10.3 , -0.02   , 3.   , -1., 0.000432, -15.43, 0.001, -7.123  }};

    int size_pos = positions_x.size();

    UNSCOPED_INFO("Test case broken!");
    REQUIRE(size_pos == positions_y.size());
    REQUIRE(size_pos == positions_z.size());

    // Define a regular grid in Galactic cartesian coordinates (units are kpc)
    const std::array<int, 3> shape {{4, 3, 2}};
    const std::array<double, 3> refpoint {{-4., 0.1, -0.3}};
    const std::array<double, 3> increment {{2.1, 0.3, 1.}};


    // Dictionaries of all models, uniform an helix are tested elsewhere
    // using pointer here since RegularField is abstract
    std::map <std::string, std::shared_ptr<RegularScalarField>> models_w_empty_constructor;
    std::map <std::string, std::shared_ptr<RegularScalarField>> models_w_regular_constructor;
    std::map <std::string, std::shared_ptr<RegularScalarField>> models_w_irregular_constructor;

    models_w_empty_constructor["YMW16"] = std::shared_ptr<YMW16> (new YMW16());
    models_w_regular_constructor["YMW16"] = std::shared_ptr<YMW16> (new YMW16(shape, refpoint, increment));
	models_w_irregular_constructor["YMW16"] = std::shared_ptr<YMW16> (new YMW16(positions_x, positions_y, positions_z));


    auto empty_model_iter = models_w_empty_constructor.begin();
    auto regular_model_iter = models_w_empty_constructor.begin();
    auto irregular_model_iter = models_w_empty_constructor.begin();

    while (empty_model_iter != models_w_empty_constructor.end()) { 
        auto pos_x_iter = positions_x.begin();
        auto pos_y_iter = positions_y.begin();
        auto pos_z_iter = positions_z.begin();
        while (pos_x_iter != positions_x.end()) { 
            UNSCOPED_INFO("RegularModels testing: at_position call failed for: ");
            CAPTURE(irregular_model_iter->first, *pos_x_iter, *pos_y_iter, *pos_z_iter);
            REQUIRE_NOTHROW((*models_w_empty_constructor[empty_model_iter->first]).at_position(*pos_x_iter, *pos_y_iter, *pos_z_iter));
            UNSCOPED_INFO("RegularModels testing: at_position contained NaN: ");
            CAPTURE(irregular_model_iter->first, *pos_x_iter, *pos_y_iter, *pos_z_iter);
            auto eval = (*models_w_empty_constructor[empty_model_iter->first]).at_position(*pos_x_iter, *pos_y_iter, *pos_z_iter);
            REQUIRE(eval != NAN);
            
            //#if autodiff_FOUND
            // test derivative
            //UNSCOPED_INFO("RegularModels testing: derivative call failed for: ");
            //CAPTURE(irregular_model_iter->first, pos_x_iter, pos_y_iter, pos_z_iter);
            //REQUIRE_NOTHROW((*models_w_empty_constructor[empty_model_iter->first]).derivative(*pos_x_iter, *pos_y_iter, *pos_z_iter));
            //#endif

            ++pos_x_iter;
            ++pos_y_iter;
            ++pos_z_iter;
            }

        UNSCOPED_INFO("RegularModels testing: on grid failed for: ");
        CAPTURE(irregular_model_iter->first);
        //std::cout << irregular_model_iter->first << std::endl;
        
        REQUIRE_THROWS_WITH((*models_w_empty_constructor[empty_model_iter->first]).on_grid(), "The class has not been initialized with a grid, hence on_grid can only be called with a grid provided.");
        double* eval_irregular = (*models_w_irregular_constructor[irregular_model_iter->first]).
        on_grid(); 
        double* eval_regular = (*models_w_regular_constructor[regular_model_iter->first]).on_grid(); 
        double* eval_no_grid_regular = (*models_w_empty_constructor[empty_model_iter->first]).on_grid(shape, refpoint, increment); 
        double* eval_no_grid_irregular = (*models_w_empty_constructor[empty_model_iter->first]).on_grid(positions_x, positions_y, positions_z);         

        REQUIRE_THAT(eval_regular, EqualsPointer(eval_no_grid_regular, 4*3*2));
        REQUIRE_THAT(eval_irregular, EqualsPointer(eval_no_grid_irregular, size_pos));
        
        ++empty_model_iter;
        ++regular_model_iter;
        ++irregular_model_iter;
    }
}

