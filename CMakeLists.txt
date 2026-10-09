cmake_minimum_required(VERSION 3.16)

project(ImagineModelsProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)

option(BUILD_PYTHON_PACKAGE "Build the Python bindings" OFF)

if(BUILD_PYTHON_PACKAGE)
    set(IMAGINE_SHARED_LIBRARY OFF CACHE BOOL "" FORCE)
endif()

add_subdirectory(${PROJECT_SOURCE_DIR}/c_library)

if(BUILD_PYTHON_PACKAGE)
    set(PYBIND11_FINDPYTHON ON)
    find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

    set(INTERNAL_PYBIND11_DIR ${PROJECT_SOURCE_DIR}/extern/pybind11)
    if(NOT EXISTS ${INTERNAL_PYBIND11_DIR}/CMakeLists.txt)
        message(FATAL_ERROR "The pybind11 submodule was not initialized! Please run 'git submodule update --init'.")
    endif()
    add_subdirectory(${INTERNAL_PYBIND11_DIR})

    set(PW ${PROJECT_SOURCE_DIR}/python_wrapper)
    set(BINDING_SOURCES
        ${PW}/module.cc
        ${PW}/grids.cc
        ${PW}/interpolation.cc
        ${PW}/regular/bases.cc
        ${PW}/regular/archimedes.cc
        ${PW}/regular/fauvet.cc
        ${PW}/regular/han.cc
        ${PW}/regular/helix.cc
        ${PW}/regular/hmr.cc
        ${PW}/regular/jaffe.cc
        ${PW}/regular/kst24.cc
        ${PW}/regular/ne2025.cc
        ${PW}/regular/plane_parallel.cc
        ${PW}/regular/yt20.cc
        ${PW}/regular/jf12.cc
        ${PW}/regular/ps.cc
        ${PW}/regular/pshirkov.cc
        ${PW}/regular/stanev.cc
        ${PW}/regular/sun.cc
        ${PW}/regular/svt22.cc
        ${PW}/regular/tf17.cc
        ${PW}/regular/tt.cc
        ${PW}/regular/uf24.cc
        ${PW}/regular/uniform.cc
        ${PW}/regular/wmap.cc
        ${PW}/regular/ymw16.cc
        ${PW}/regular/xh24.cc
    )
    if(IMAGINE_HAS_FFTW)
        list(APPEND BINDING_SOURCES
            ${PW}/random/bases.cc
            ${PW}/random/es_random.cc
            ${PW}/random/uf26_random.cc
            ${PW}/random/sun_random.cc
            ${PW}/random/jaffe_random.cc
            ${PW}/random/orlando26_random.cc
            ${PW}/random/gaussian_scalar.cc
            ${PW}/random/jf12_random.cc
            ${PW}/random/log_normal.cc
        )
    endif()

    pybind11_add_module(_core ${BINDING_SOURCES})
    target_link_libraries(_core PRIVATE ImagineModels::ImagineModels)
    set_target_properties(_core PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

    install(TARGETS _core LIBRARY DESTINATION ImagineModels COMPONENT python)
endif()
