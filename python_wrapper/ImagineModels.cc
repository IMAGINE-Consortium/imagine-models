#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace pybind11::literals;


//PYBIND11_MAKE_OPAQUE(std::vector<double>);
//PYBIND11_MAKE_OPAQUE(std::array<double, 3>);
//PYBIND11_MAKE_OPAQUE(std::array<int, 3>);

#if autodiff_FOUND
  #include "include/autodiff_wrapper.h"
#endif

#include "include/FieldBases.h"

#include "include/regular/RegularFieldBases.h"
#include "include/regular/SunWrapper.h"
#include "include/regular/HanWrapper.h"
#include "include/regular/WMAPWrapper.h"
#include "include/regular/StanevBSSWrapper.h"
#include "include/regular/FauvetWrapper.h"
#include "include/regular/TinyakovTkachevWrapper.h"
#include "include/regular/HarariMollerachRouletWrapper.h"
#include "include/regular/UniformWrapper.h"
#include "include/regular/YMW16Wrapper.h"
#include "include/regular/HelixWrapper.h"
#include "include/regular/JaffeWrapper.h"
#include "include/regular/HarariMollerachRouletWrapper.h"
#include "include/regular/TinyakovTkachevWrapper.h"
#include "include/regular/RegularJF12Wrapper.h"

#if FFTW_FOUND
  #include "include/random/RandomFieldBases.h"
  #include "include/random/RandomJF12Wrapper.h"
  #include "include/random/EnsslinSteiningerWrapper.h"
  #include "include/random/GaussianScalarWrapper.h"
  #include "include/random/LogNormalWrapper.h"
#endif

void FieldBases(py::module_ &);
void RegularFieldBases(py::module_ &);

void RegularJF12(py::module_ &);
void Helix(py::module_ &);
void Uniform(py::module_ &);
void Jaffe(py::module_ &);
void Sun2008(py::module_ &);
void Han2018(py::module_ &);
void StanevBSS(py::module_ &);
void TinyakovTkachev(py::module_ &);
void HarariMollerachRoulet(py::module_ &);
void WMAP(py::module_ &);
void Fauvet(py::module_ &);

void YMW(py::module_ &);

#if FFTW_FOUND
    void RandomFieldBases(py::module_ &);
    
    void RandomJF12(py::module_ &);
    void EnsslinSteininger(py::module_ &);
    void GaussianScalar(py::module_ &);
    void LogNormal(py::module_ &);
#endif


PYBIND11_MODULE(_ImagineModels, m) {
    m.doc() = "IMAGINE Model Library";

    FieldBases(m);
    RegularFieldBases(m);
    RegularJF12(m);
    Jaffe(m);
    Uniform(m);
    Sun2008(m);
    Han2018(m);
    StanevBSS(m);
    Helix(m);
    TinyakovTkachev(m);
    HarariMollerachRoulet(m);
    YMW(m);
    WMAP(m);
    Fauvet(m);

    #if FFTW_FOUND
      RandomFieldBases(m);    
      RandomJF12(m);
      EnsslinSteininger(m);
      GaussianScalar(m);
      LogNormal(m);
    #endif

}