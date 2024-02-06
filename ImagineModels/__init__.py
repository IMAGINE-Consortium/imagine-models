from _ImagineModels import VectorFieldBase, ScalarFieldBase, RegularVectorField, RegularScalarField, \
    JF12RegularField,  JaffeMagneticField, HelixMagneticField, UniformMagneticField, UniformDensityField, YMW16, SunMagneticField, \
    HanMagneticField,\
    WMAPMagneticField, TTMagneticField, HMRMagneticField, FauvetMagneticField, StanevBSSMagneticField, TFMagneticField, PshirkovMagneticField, ArchimedeanMagneticField

    
try: 
    from _ImagineModels import JF12RandomField, ESRandomField, GaussianScalarField, LogNormalScalarField
    __has_random_fields__ = True
except ImportError:
    __has_random_fields__ = False
    
try:  
    _jf12 = JF12RegularField()
    _ = _jf12.active_diff
    __has_autodiff__ = True
except AttributeError:
    __has_autodiff__ = False
    
    
        
from .HelperFunctions.CoordinateConversions import cyl2cart

from .MagneticFields.RegularMagneticFields import AxiSymmetricSpiral
