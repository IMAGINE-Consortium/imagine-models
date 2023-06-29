#include <functional>
#include <cmath>

#include "ImagineModels/Field.h"
#include "ImagineModels/RegularField.h"


//Sun et al. A&A V.477 2008 ASS+RING model magnetic field
class Sun2008MagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;
        double b_B0 = 2.;
        double b_Rsun = 8.5; 
        double b_R0 = 10.;
        double b_z0 = 1.;
        double b_Rc = 5.;
        double b_Bc = 2.;
        double b_pitch_deg = -12.;
        
        double bH_B0 = 2.;
        double bH_z0 = 1.5;
        double bH_z1a = 0.2;
        double bH_z1b = 0.4;
        double bH_R0 = 4.;

        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const;
 };
