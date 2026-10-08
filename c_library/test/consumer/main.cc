#include <cmath>
#include <iostream>

#include "ImagineModels/ImagineModels.h"

int main() {
    imagine::JF12MagneticField jf12;
    imagine::Vec3<double> b = jf12.at_position(-8.5, 0., 0.1);
    std::cout << "ImagineModels " << IMAGINE_VERSION << ": JF12 at sun " << double(b[0]) << " " << double(b[1]) << " "
              << double(b[2]) << std::endl;
    if (!std::isfinite(double(b[0])))
        return 1;
#if IMAGINE_HAS_AUTODIFF
    std::cout << "derivative columns: " << jf12.derivative(-8.5, 0., 0.1).cols() << std::endl;
#endif
#if IMAGINE_HAS_FFTW
    imagine::GaussianScalarField gauss;
    imagine::ScalarGridData grid = gauss.sample(imagine::RegularGrid({4, 4, 4}, {0., 0., 0.}, {1., 1., 1.}), 3);
    std::cout << "random grid value: " << grid(0, 0) << std::endl;
#endif
    return 0;
}
