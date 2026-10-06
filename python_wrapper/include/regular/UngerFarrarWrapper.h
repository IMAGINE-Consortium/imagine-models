#ifndef UNGERFARRARWRAPPER_H
#define UNGERFARRARWRAPPER_H

#include "ImagineModels/UngerFarrar.h"
#include "model_bindings.h"

void UF24(py::module_ &m)
{
    bind_regular_model<UFMagneticField>(m, "UFMagneticField")
        .def_readwrite("activeModel", &UFMagneticField::activeModel)
        .def_readonly("possibleModels", &UFMagneticField::possibleModels)
        .def_readonly("all_parameters", &UFMagneticField::all_parameters)
        .def("set_parameters", [](UFMagneticField &self, std::string model_type) {
            self.set_parameters(model_type); 
        });
}

#endif
