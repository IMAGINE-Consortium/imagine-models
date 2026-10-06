#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/stl_bind.h>

#include "ImagineModels/ImagineModels.h"

#if IMAGINE_HAS_AUTODIFF
  #include <pybind11/eigen.h>
  #include "include/autodiff_wrapper.h"
#endif

namespace py = pybind11;
using namespace pybind11::literals;
using namespace imagine;

#include "include/grid_bindings.h"

#include "include/regular/RegularFieldBases.h"
#include "include/regular/SunWrapper.h"
#include "include/regular/HanWrapper.h"
#include "include/regular/PshirkovWrapper.h"
#include "include/regular/ArchimedeanWrapper.h"
#include "include/regular/WMAPWrapper.h"
#include "include/regular/StanevBSSWrapper.h"
#include "include/regular/FauvetWrapper.h"
#include "include/regular/TinyakovTkachevWrapper.h"
#include "include/regular/HarariMollerachRouletWrapper.h"
#include "include/regular/UniformWrapper.h"
#include "include/regular/YMW16Wrapper.h"
#include "include/regular/TF17Wrapper.h"
#include "include/regular/HelixWrapper.h"
#include "include/regular/JaffeWrapper.h"
#include "include/regular/HarariMollerachRouletWrapper.h"
#include "include/regular/TinyakovTkachevWrapper.h"
#include "include/regular/RegularJF12Wrapper.h"
#include "include/regular/UngerFarrarWrapper.h"
#include "include/regular/SVT22Wrapper.h"


#if IMAGINE_HAS_FFTW
#include "include/random/RandomFieldBases.h"
#include "include/random/RandomJF12Wrapper.h"
#include "include/random/EnsslinSteiningerWrapper.h"
#include "include/random/GaussianScalarWrapper.h"
#include "include/random/LogNormalWrapper.h"
#endif

void Grids(py::module_ &);
void RegularFieldBases(py::module_ &);

void RegularJF12(py::module_ &);
void Helix(py::module_ &);
void Uniform(py::module_ &);
void Jaffe(py::module_ &);
void Archimedes(py::module_ &);
void Pshirkov(py::module_ &);
void Sun2008(py::module_ &);
void TF17(py::module_ &);
void UF24(py::module_ &);
void Han2018(py::module_ &);
void StanevBSS(py::module_ &);
void TinyakovTkachev(py::module_ &);
void HarariMollerachRoulet(py::module_ &);
void WMAP(py::module_ &);
void Fauvet(py::module_ &);
void YMW(py::module_ &);
void SVT22(py::module_ &);

#if IMAGINE_HAS_FFTW
void RandomFieldBases(py::module_ &);

void RandomJF12(py::module_ &);
void EnsslinSteininger(py::module_ &);
void GaussianScalar(py::module_ &);
void LogNormal(py::module_ &);

#endif

PYBIND11_MODULE(_ImagineModels, m)
{
  m.doc() = "IMAGINE Model Library";

  Grids(m);
  RegularFieldBases(m);
  RegularJF12(m);
  Jaffe(m);
  Uniform(m);
  Sun2008(m);
  Han2018(m);
  TF17(m);
  UF24(m);
  Archimedes(m);
  Pshirkov(m);
  StanevBSS(m);
  Helix(m);
  TinyakovTkachev(m);
  HarariMollerachRoulet(m);
  YMW(m);
  WMAP(m);
  Fauvet(m);
  SVT22(m);
#if IMAGINE_HAS_FFTW
  RandomFieldBases(m);
  RandomJF12(m);
  EnsslinSteininger(m);
  GaussianScalar(m);
  LogNormal(m);
#endif
}