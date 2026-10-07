from ._core import __version__, has_autodiff, has_fftw
from ._core import RegularGrid, IrregularGrid, PointCloud, GridError, RegularVectorField, RegularScalarField, SVT22, JF12RegularField,  JaffeMagneticField, HelixMagneticField, UniformMagneticField, UniformDensityField, YMW16, SunMagneticField, HanMagneticField,WMAPMagneticField, TTMagneticField, HMRMagneticField, FauvetMagneticField, StanevBSSMagneticField, TFMagneticField, PshirkovMagneticField, ArchimedeanMagneticField, UFMagneticField, XH24MagneticField

if has_fftw:
    from ._core import RandomVectorField, RandomScalarField, JF12RandomField, ESRandomField, UF26RandomField, GaussianScalarField, LogNormalScalarField

__has_random_fields__ = has_fftw
__has_autodiff__ = has_autodiff

from .HelperFunctions.CoordinateConversions import cyl2cart, cart2cyl
from .HelperFunctions.PlottingHelpers import plot_slice

from .MagneticFields.RegularMagneticFields import AxiSymmetricSpiral
