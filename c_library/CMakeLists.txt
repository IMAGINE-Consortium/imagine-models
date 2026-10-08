cmake_minimum_required(VERSION 3.16)

project(ImagineModels VERSION 0.2 DESCRIPTION "A library for Galactic models." LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    include(CTest)
endif()

option(USE_AUTODIFF "Enable derivatives w.r.t. model parameters (requires autodiff and Eigen3)" ON)
option(USE_FFTW "Enable random field models (requires FFTW3)" ON)
option(IMAGINE_SHARED_LIBRARY "Build ImagineModels as a shared library" ON)

set(IMAGINE_HAS_AUTODIFF 0)
if(USE_AUTODIFF)
    find_package(autodiff QUIET)
    find_package(Eigen3 QUIET NO_MODULE)
    if(autodiff_FOUND AND Eigen3_FOUND)
        set(IMAGINE_HAS_AUTODIFF 1)
        message(STATUS "ImagineModels: autodiff and Eigen3 found, derivatives enabled.")
    else()
        message(STATUS "ImagineModels: autodiff or Eigen3 not found, derivatives disabled.")
    endif()
else()
    message(STATUS "ImagineModels: derivatives disabled (USE_AUTODIFF=OFF).")
endif()

set(IMAGINE_HAS_FFTW 0)
set(FINDFFTW_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../extern/findFFTW)
if(USE_FFTW)
    if(NOT EXISTS ${FINDFFTW_DIR}/FindFFTW.cmake)
        message(FATAL_ERROR "The findFFTW submodule was not initialized! Please run 'git submodule update --init'.")
    endif()
    list(APPEND CMAKE_MODULE_PATH ${FINDFFTW_DIR})
    find_package(FFTW QUIET COMPONENTS DOUBLE_LIB)
    if(FFTW_FOUND)
        set(IMAGINE_HAS_FFTW 1)
        message(STATUS "ImagineModels: FFTW found (${FFTW_DOUBLE_LIB}), random fields enabled.")
    else()
        message(STATUS "ImagineModels: FFTW not found, random fields disabled.")
    endif()
else()
    message(STATUS "ImagineModels: random fields disabled (USE_FFTW=OFF).")
endif()

set(IM_HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/headers)
set(IM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source)
set(IM_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)

set(IMAGINE_HAS_AUTODIFF ${IMAGINE_HAS_AUTODIFF} CACHE INTERNAL "")
set(IMAGINE_HAS_FFTW ${IMAGINE_HAS_FFTW} CACHE INTERNAL "")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${IM_GENERATED_DIR}/ImagineModels/config.h @ONLY)

set(IM_SRC_FILES
    ${IM_SOURCE_DIR}/helix.cc
    ${IM_SOURCE_DIR}/jaffe.cc
    ${IM_SOURCE_DIR}/kst24.cc
    ${IM_SOURCE_DIR}/ne2025.cc
    ${IM_SOURCE_DIR}/ne2025_tables.cc
    ${IM_SOURCE_DIR}/plane_parallel.cc
    ${IM_SOURCE_DIR}/yt20.cc
    ${IM_SOURCE_DIR}/jf12.cc
    ${IM_SOURCE_DIR}/ymw16.cc
    ${IM_SOURCE_DIR}/xh24.cc
    ${IM_SOURCE_DIR}/tf17.cc
    ${IM_SOURCE_DIR}/sun.cc
    ${IM_SOURCE_DIR}/pshirkov.cc
    ${IM_SOURCE_DIR}/archimedes.cc
    ${IM_SOURCE_DIR}/han.cc
    ${IM_SOURCE_DIR}/wmap.cc
    ${IM_SOURCE_DIR}/stanev.cc
    ${IM_SOURCE_DIR}/fauvet.cc
    ${IM_SOURCE_DIR}/tt.cc
    ${IM_SOURCE_DIR}/hmr.cc
    ${IM_SOURCE_DIR}/svt22.cc
    ${IM_SOURCE_DIR}/uf24.cc
)

if(IMAGINE_HAS_FFTW)
    list(APPEND IM_SRC_FILES
        ${IM_SOURCE_DIR}/random_field.cc
        ${IM_SOURCE_DIR}/random_scalar_field.cc
        ${IM_SOURCE_DIR}/random_vector_field.cc
        ${IM_SOURCE_DIR}/gaussian_scalar.cc
        ${IM_SOURCE_DIR}/log_normal.cc
        ${IM_SOURCE_DIR}/jf12_random.cc
        ${IM_SOURCE_DIR}/es_random.cc
        ${IM_SOURCE_DIR}/uf26_random.cc
    )
endif()

if(IMAGINE_SHARED_LIBRARY)
    add_library(ImagineModels SHARED ${IM_SRC_FILES})
else()
    add_library(ImagineModels STATIC ${IM_SRC_FILES})
    set_target_properties(ImagineModels PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
add_library(ImagineModels::ImagineModels ALIAS ImagineModels)

target_compile_features(ImagineModels PUBLIC cxx_std_17)
target_include_directories(ImagineModels PUBLIC
    $<BUILD_INTERFACE:${IM_HEADER_DIR}>
    $<BUILD_INTERFACE:${IM_GENERATED_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if(IMAGINE_HAS_AUTODIFF)
    target_link_libraries(ImagineModels PUBLIC autodiff::autodiff Eigen3::Eigen)
endif()
if(IMAGINE_HAS_FFTW)
    target_link_libraries(ImagineModels PUBLIC FFTW::Double)
endif()

set_target_properties(ImagineModels PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)

if(BUILD_TESTING)
    find_package(Catch2 3 QUIET)
    if(NOT Catch2_FOUND)
        include(FetchContent)
        FetchContent_Declare(Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG v3.5.2
            GIT_SHALLOW TRUE)
        FetchContent_MakeAvailable(Catch2)
        list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    endif()
    include(Catch)

    set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)
    set(TEST_SOURCES
        ${TEST_DIR}/regular/test_grid.cc
        ${TEST_DIR}/regular/test_interpolation.cc
        ${TEST_DIR}/regular/test_parameters.cc
        ${TEST_DIR}/regular/test_positions.cc
        ${TEST_DIR}/regular/test_models.cc
        ${TEST_DIR}/regular/test_derivatives.cc
    )
    if(IMAGINE_HAS_FFTW)
        list(APPEND TEST_SOURCES ${TEST_DIR}/random/test_random.cc)
    endif()
    add_executable(imagine_tests ${TEST_SOURCES})
    target_include_directories(imagine_tests PRIVATE ${TEST_DIR})
    target_link_libraries(imagine_tests PRIVATE ImagineModels Catch2::Catch2WithMain)
    catch_discover_tests(imagine_tests PROPERTIES SKIP_RETURN_CODE 4)
endif()

set(IM_CMAKE_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/ImagineModels)

install(TARGETS ImagineModels EXPORT ImagineModelsTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY ${IM_HEADER_DIR}/ImagineModels DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${IM_GENERATED_DIR}/ImagineModels/config.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ImagineModels)
if(IMAGINE_HAS_FFTW)
    install(DIRECTORY ${IM_HEADER_DIR}/ImagineModelsRandom DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(FILES ${FINDFFTW_DIR}/FindFFTW.cmake DESTINATION ${IM_CMAKE_INSTALL_DIR})
endif()

install(EXPORT ImagineModelsTargets NAMESPACE ImagineModels:: DESTINATION ${IM_CMAKE_INSTALL_DIR})
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/ImagineModelsConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ImagineModelsConfig.cmake
    INSTALL_DESTINATION ${IM_CMAKE_INSTALL_DIR}
)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/ImagineModelsConfigVersion.cmake
    COMPATIBILITY SameMinorVersion
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/ImagineModelsConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ImagineModelsConfigVersion.cmake
    DESTINATION ${IM_CMAKE_INSTALL_DIR}
)

set(IM_PC_REQUIRES "")
if(IMAGINE_HAS_AUTODIFF)
    string(APPEND IM_PC_REQUIRES " eigen3")
endif()
if(IMAGINE_HAS_FFTW)
    string(APPEND IM_PC_REQUIRES " fftw3")
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ImagineModels.pc.in ${CMAKE_CURRENT_BINARY_DIR}/ImagineModels.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ImagineModels.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
