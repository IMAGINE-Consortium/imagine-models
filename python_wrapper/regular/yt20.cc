#include "ImagineModels/YT20.h"
#include "../bindings.h"
#include "../model_bindings.h"

void bind_yt20(py::module_ &m) {
    bind_regular_model<YT20>(m, "YT20").def_readwrite("mu", &YT20::mu).def_readwrite("mu_e", &YT20::mu_e);
}
