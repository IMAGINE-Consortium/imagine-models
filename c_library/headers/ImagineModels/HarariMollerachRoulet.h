#include <functional>
#include <cmath>

#include "ImagineModels/Field.h"
#include "ImagineModels/RegularField.h"


// Harari, Mollerach, Roulet (HMR) see https://arxiv.org/abs/astro-ph/9906309, implementation of https://arxiv.org/pdf/astro-ph/0510444.pdf 

class HMRMagneticField : public RegularVectorField  {
    protected:
        bool DEBUG = false;
    public:
        using RegularVectorField :: RegularVectorField;
        double b_Rsun = 8.5; // kpc
        double b_r_max = 20.; // kpc
        double b_z1 = 0.3; // kpc
        double b_z2 = 4.; // kpc
        double b_r1 = 2.; // kpc
        double b_p = -10; // degree
        double b_epsilon0 = 10.55; // kpc


        std::array<double, 3>  at_position (const double &x, const double &y, const double &z) const;
 };