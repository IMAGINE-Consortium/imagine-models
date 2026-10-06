cmake_minimum_required(VERSION 3.16)

project(_ImagineModels DESCRIPTION "A Python module for Galactic models.")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

include(CTest)
include(GNUInstallDirs)

option(BUILD_PYTHON_PACKAGE "Build the Python bindings" OFF)

### INSTALLATION 

## C++

# now we add the c_library, the CMakeLists in there will install the ImagineModels library 
add_subdirectory(${PROJECT_SOURCE_DIR}/c_library)

## Python, if requested

if(BUILD_PYTHON_PACKAGE)

    # Checking for Python 3
    set(INTERNAL_PYBIND11_DIR "${PROJECT_SOURCE_DIR}/extern/pybind11")
    set(Python_FIND_VIRTUALENV FIRST)
    find_package(Python3 COMPONENTS Interpreter Development)
    if(Python3_FOUND)
        message("-- Python ${Python3_VERSION} found. Enabling support for pybind.")
        message("Python executable: ${Python3_EXECUTABLE}")

        ## looking for pybind11
        if(NOT EXISTS ${INTERNAL_PYBIND11_DIR})
            message(FATAL_ERROR "The pybind11 submodule was not initialized! Please run 'git submodule update --init'." )
        endif()
        add_subdirectory(${INTERNAL_PYBIND11_DIR})

        ## Here, we add the main binding code
        pybind11_add_module(_ImagineModels ${PROJECT_SOURCE_DIR}/python_wrapper/ImagineModels.cc)

        set_target_properties(_ImagineModels PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
        set_target_properties(_ImagineModels PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

        target_link_libraries(_ImagineModels PRIVATE ImagineModels::ImagineModels)

        if(UNIX AND NOT APPLE)
            message("-- _ImagineModels: We are on a non MacOS unix system, installing to ${Python3_SITEARCH}") 
            set_target_properties(_ImagineModels PROPERTIES INSTALL_RPATH "${Python3_SITEARCH}")
        elseif(APPLE)
        message("-- _ImagineModels: We are on a MacOS system, installing to loader_path")
            set_target_properties(_ImagineModels PROPERTIES INSTALL_RPATH "@loader_path")
        else() 
            message(FATAL_ERROR "Only unix based systems are supported")
        endif()
        
        # EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
        # define (VERSION_INFO) here.
        target_compile_definitions(_ImagineModels PRIVATE VERSION_INFO=${PROJECT_VERSION})
        set(install_path "${Python3_SITEARCH}")

        
    else()
        message("-- Could not find the Python Interpreter! If you want to build the python package, please provide provide a hint to CMake, using '-D CMAKE_PYTHON_EXECUTABLE'.")
    endif()
endif()