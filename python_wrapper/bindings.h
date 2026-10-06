#ifndef IMAGINE_BINDINGS_H
#define IMAGINE_BINDINGS_H

#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "ImagineModels/config.h"
#include "ImagineModels/Grid.h"

#if IMAGINE_HAS_AUTODIFF
#include <pybind11/eigen.h>
#endif

namespace py = pybind11;
using namespace pybind11::literals;
using namespace imagine;

template <int N>
py::array_t<double> to_numpy(GridData<N> &&grid_data) {
  auto owned = new GridData<N>(std::move(grid_data));
  py::capsule owner(owned, [](void *p) { delete static_cast<GridData<N> *>(p); });
  std::vector<py::ssize_t> shape;
  if (N > 1)
    shape.push_back(N);
  for (int s : owned->shape)
    shape.push_back(s);
  return py::array_t<double>(shape, owned->data.data(), owner);
}

void bind_grids(py::module_ &m);
void bind_regular_bases(py::module_ &m);
void bind_archimedes(py::module_ &m);
void bind_fauvet(py::module_ &m);
void bind_han(py::module_ &m);
void bind_hmr(py::module_ &m);
void bind_helix(py::module_ &m);
void bind_jaffe(py::module_ &m);
void bind_pshirkov(py::module_ &m);
void bind_jf12(py::module_ &m);
void bind_stanev(py::module_ &m);
void bind_sun(py::module_ &m);
void bind_svt22(py::module_ &m);
void bind_tf17(py::module_ &m);
void bind_tt(py::module_ &m);
void bind_uf24(py::module_ &m);
void bind_uniform(py::module_ &m);
void bind_wmap(py::module_ &m);
void bind_ymw16(py::module_ &m);

#if IMAGINE_HAS_FFTW
void bind_random_bases(py::module_ &m);
void bind_es_random(py::module_ &m);
void bind_gaussian_scalar(py::module_ &m);
void bind_lognormal(py::module_ &m);
void bind_jf12_random(py::module_ &m);
#endif

#endif
