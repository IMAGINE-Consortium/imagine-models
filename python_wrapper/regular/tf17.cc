#include "../bindings.h"
#include "ImagineModels/TF17.h"
#include "../model_bindings.h"

void bind_tf17(py::module_ &m)
{
    bind_regular_model<TFMagneticField>(m, "TFMagneticField")
        .def_readwrite("activeDiskModel", &TFMagneticField::activeDiskModel)
        .def_readwrite("activeHaloModel", &TFMagneticField::activeHaloModel)
        .def_readonly("possibleDiskModels", &TFMagneticField::possibleDiskModels)
        .def_readonly("possibleHaloModels", &TFMagneticField::possibleHaloModels)
        .def("set_params", [](TFMagneticField &self, std::string dtype, std::string htype) {
            self.set_params(dtype, htype); 
        });
}
