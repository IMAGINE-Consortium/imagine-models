from ._core import (
    YMW16,
    ArchimedeanMagneticField,
    FauvetMagneticField,
    GridError,
    HanMagneticField,
    HelixMagneticField,
    HMRMagneticField,
    IrregularGrid,
    JaffeMagneticField,
    JF12MagneticField,
    KST24MagneticField,
    PointCloud,
    PshirkovMagneticField,
    RegularGrid,
    RegularScalarField,
    RegularVectorField,
    StanevBSSMagneticField,
    SunMagneticField,
    SVT22MagneticField,
    TF17MagneticField,
    TTMagneticField,
    UF24MagneticField,
    UniformDensityField,
    UniformMagneticField,
    WMAPMagneticField,
    XH24MagneticField,
    __version__,
    has_autodiff,
    has_fftw,
    interpolate,
)

if has_fftw:
    from ._core import (
        ESRandomField,
        GaussianScalarField,
        JF12RandomField,
        LogNormalScalarField,
        RandomScalarField,
        RandomVectorField,
        UF26RandomField,
    )

__has_random_fields__ = has_fftw
__has_autodiff__ = has_autodiff

from .HelperFunctions.CoordinateConversions import cart2cyl, cyl2cart
from .HelperFunctions.PlottingHelpers import plot_slice
from .MagneticFields.RegularMagneticFields import AxiSymmetricSpiral
