#ifndef ARRAY_CONVERTERS_H
#define ARRAY_CONVERTERS_H


namespace py = pybind11;
using namespace pybind11::literals;


// helper functions to avoid making a copy when returning a py::array_t
// based on
// author: https://github.com/YannickJadoul
// source: https://github.com/pybind/pybind11/issues/1042#issuecomment-642215028
template <typename Sequence>
inline py::array_t<typename Sequence::value_type> as_pyarray(Sequence &&seq) {
  auto size = seq.size();
  auto data = seq.data();
  std::unique_ptr<Sequence> seq_ptr =
      std::make_unique<Sequence>(std::move(seq));
  auto capsule = py::capsule(seq_ptr.get(), [](void *p) {
    std::unique_ptr<Sequence>(reinterpret_cast<Sequence *>(p));
  });
  seq_ptr.release();
  return py::array(size, data, capsule);
}

inline py::array_t<double> from_pointer_to_pyarray(double* data, size_t arr_size_x, size_t arr_size_y, size_t arr_size_z) {
  
  py::capsule capsule(data, [](void *f) {
      std::unique_ptr<double>(reinterpret_cast<double*>(f));
      });

  size_t arr_size = arr_size_x*arr_size_y*arr_size_z;
  py::array_t<double> arr = py::array(arr_size, data, capsule);

  return arr.reshape({arr_size_x, arr_size_y, arr_size_z});
}


inline py::list from_pointer_array_to_list_pyarray(std::array<double*, 3> seq, size_t arr_size_x, size_t arr_size_y, size_t arr_size_z) {

  py::list li;

  for (int i = 0; i<3; ++i) {
    py::array_t<double> arr = from_pointer_to_pyarray(std::move(seq[i]), arr_size_x, arr_size_y, arr_size_z);

    li.append(arr);
  }
  return li;
}

#endif