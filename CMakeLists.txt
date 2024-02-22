cmake_minimum_required(VERSION 3.16)

project(_ImagineModels DESCRIPTION "A Python module for Galactic models.")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

include(CTest)
include(FetchContent)
include(GNUInstallDirs)

set(USE_AUTODIFF ON CACHE BOOL OFF)
set(COMPILE_TESTS ON CACHE BOOL OFF)
set(USE_FFTW ON CACHE BOOL OFF)
set(BUILD_PYTHON_PACKAGE OFF CACHE BOOL OFF)

### DEPENDENCIES
# First we try to find dependecies, if requested.

## starting with the tests

if (COMPILE_TESTS)
    message("-- Trying to find Catch2 to compile the test cases. Turn off via setting 'COMPILE_TESTS=OFF' in your environment before installation. ")
    find_package(Catch2 3)
    if (Catch2_FOUND)         
        message("-- global installation of catch2 library found. Tests will be compiled.")
    else()
    message("-- no global installation of catch2 library found. Trying to include it as a submodule.")
        set(INTERNAL_findCatch2_DIR ${PROJECT_SOURCE_DIR}/extern/Catch2)
        if(NOT EXISTS ${INTERNAL_findCatch2_DIR})
            message(FATAL_ERROR "The catch2 submodule was not initialized! Please run 'git submodule update --init'. Alternativly you can install catch2 globally on your system." )
        endif()
        add_subdirectory(${INTERNAL_findCatch2_DIR})
        message("-- catch2 library found as a submodule. Tests will be compiled.")
    endif()
else()
    set(Catch2_FOUND OFF) 
    message("-- testing manually disabled")
endif()

## then AUTODIFF

if (USE_AUTODIFF)
    message("-- Trying to find autodiff. Turn off via setting 'USE_AUTODIFF=OFF' in your environment before installation. ")
    find_package(autodiff)
    if (autodiff_FOUND)         
        message("-- autodiff library found: attempting to find Eigen3.")
        set(autodiff_INCLUDE_DIRS /usr/include/eigen3)
    
        find_package(Eigen3 REQUIRED NO_MODULE)

        if(NOT Eigen3_FOUND)
            message("-- Could not find the Eigen3 library! Autodifferention disabled.")
        else()
            message("-- Eigen3 library found: Autodifferention enabled.")
    endif()
    else()
        message("-- Could not find the autodiff library! Autodifferention disabled.")
    endif()
else()
    set(autodiff_FOUND OFF) 
    message("-- autodiff library manually disabled")
endif()

add_compile_definitions(autodiff_FOUND=${autodiff_FOUND})
                
## and FFTW

if (USE_FFTW)
    message("-- Trying to find FFTW for random field generation. Turn off via setting 'USE_FFTW=OFF' in your environment before installation.")
    # Activating the findFFTW cmake package
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/c_library/findFFTW")
    message("-- Trying to find FFTW. Turn off via 'cmake ..DUSE_FFTW=OFF'")
    set(INTERNAL_findFFTW_DIR ${PROJECT_SOURCE_DIR}/extern/findFFTW)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${INTERNAL_findFFTW_DIR}")
    if(NOT EXISTS ${INTERNAL_findFFTW_DIR})
        message(FATAL_ERROR "The findFFTW submodule was not initialized! Please run 'git submodule update --init'." )
    endif()
    find_package(FFTW QUIET COMPONENTS DOUBLE_LIB)
    if(FFTW_FOUND)
        message("-- FFTW libraries found: ${FFTW_DOUBLE_LIB}; will include random field generation routines.")
            include_directories(${FFTW_INCLUDE_DIRS})
            set(LIBRARIES ${LIBRARIES} ${FFTW_DOUBLE_LIB})
    else()
        message("-- FFTW libraries not found: random field generation disabled.")
    endif()
else()
    set(FFTW_FOUND OFF) 
    message("-- fftw library manually disabled")
endif()

message("-- libraries ${LIBRARIES}")

add_compile_definitions(FFTW_FOUND=${FFTW_FOUND})

### INSTALLATION 

## C++

# now we add the c_library, the CMakeLists in there will install the ImagineModels library 
add_subdirectory(${PROJECT_SOURCE_DIR}/c_library)

## Python, if requested

if(BUILD_PYTHON_PACKAGE)

    set(IM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/c_library/headers)

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

        target_include_directories(_ImagineModels PUBLIC ${IM_INCLUDE_DIR}/ImagineModels)

        ## optional dependencies
        if(FFTW_FOUND)
            target_include_directories(_ImagineModels PUBLIC ${IM_INCLUDE_DIR}/ImagineModelsRandom)
            target_include_directories(_ImagineModels PUBLIC ${FFTW_INCLUDE_DIRS})
        endif()

        set_target_properties(_ImagineModels PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
        
        set_target_properties(_ImagineModels PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

        target_link_libraries(_ImagineModels PUBLIC ImagineModels ${LIBRARIES})
        if (autodiff_FOUND)
            target_link_libraries(_ImagineModels PUBLIC autodiff::autodiff)
            target_link_libraries(_ImagineModels PUBLIC Eigen3::Eigen)
        endif()

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