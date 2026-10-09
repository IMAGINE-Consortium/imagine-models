#include "bindings.h"

PYBIND11_MODULE(_core, m) {
    m.doc() = "IMAGINE Model Library";
    m.attr("__version__") = IMAGINE_VERSION;
    m.attr("has_autodiff") = bool(IMAGINE_HAS_AUTODIFF);
    m.attr("has_fftw") = bool(IMAGINE_HAS_FFTW);

    bind_grids(m);
    bind_interpolation(m);
    bind_regular_bases(m);
    bind_archimedes(m);
    bind_fauvet(m);
    bind_han(m);
    bind_hmr(m);
    bind_helix(m);
    bind_jaffe(m);
    bind_kst24(m);
    bind_ne2025(m);
    bind_plane_parallel(m);
    bind_yt20(m);
    bind_ps(m);
    bind_pshirkov(m);
    bind_jf12(m);
    bind_stanev(m);
    bind_sun(m);
    bind_svt22(m);
    bind_tf17(m);
    bind_tt(m);
    bind_uf24(m);
    bind_uniform(m);
    bind_wmap(m);
    bind_ymw16(m);
    bind_xh24(m);
#if IMAGINE_HAS_FFTW
    bind_random_bases(m);
    bind_es_random(m);
    bind_gaussian_scalar(m);
    bind_lognormal(m);
    bind_jf12_random(m);
    bind_uf26_random(m);
    bind_sun_random(m);
    bind_jaffe_random(m);
    bind_orlando26_random(m);
#endif
}
