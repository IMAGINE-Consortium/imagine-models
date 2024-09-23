import numpy as np

from ImagineModels import RegularVectorField, cyl2cart, cart2cyl


class AxiSymmetricSpiral(RegularVectorField):
    """
            FUNCTION
                ===  ASS : AxiSymmetric Spiral Galactic magnetic field model ===
            ASS(coordinates,**kwargs{numerous parameters})
                " vec{B} = B_rho vec{u_rho} + B_phi vec{u_phi} + B_z vec{u_z}
                    where
                   " B_rho = B_amp(rho,z) * sin(pitch) * cos(Xi(z))
                   " B_phi = B_amp(rho,z) * cos(pitch) * cos(Xi(z))
                   " B_z   = B_amp(rho,z) * sin(Xi(z))
                   with
                       " Xi(z) = Xi_0 * tanh(z / z_0)
                       " B_amp(rho,z) = cst,
                                        funct of sph. radial coord
                                        funct of cyl. radial coord
                       " pitch = pitch angle
            INPUT:
            ------
             - coord :  an (3,n)-array with Galactic coordinates at which the
                        Bfield has to be evaluated. Default format is cylindrical.
             **kwargs :
               - coord_format : string that specifies the coord format of the INPUT
                                coordinates. Can be cartesian, spherical or cylindrical
                                (the output are in cylindrical, what so ever)
                                Default is 'cartesian'
               - B_amp_type : string that specifies if the amplitude is to be
                              a constante ['cst'], a function of the cyl. radial
                              coordinate ['cyl'] or a function of the SPH. radial
                              coordinate ['sph']
                              Default = 'cyl'
              [model parameter]
               - B_0 : amplitude of the regular large-scale B field at the Sun [muG]
               - B_amp_param : additional parameter for the radial dependence
                               (if it is the case)
               - B_amp_type : string to specify the radial dependence
                              'cst', 'cyl','sph' for constant, cylindircal or spherical

               - pitch : pitch angle [rad] of the spiral arms [Default = 11.5 deg]
               - rho_0 : radial scale = dist. between the Sun and the Gal Centre [kpc]
               - Xi_0 :  tilt parameter in [rad]
               - z_0 :   a vertical scale in [kpc]
               > getBFieldDefault('ASS') for default setting values of parameters

            OUTPUT:
            ------
            B_rho, B_z, B_phi : the component of the vectorial B field in CYLINDRICAL
                                coordinate system centred on the Galactic Centre at
                                all locations specified in 'coord'
            Created on Jul 4 2017
            @author: V.Pelgrims
            """

    def __init__(self, b_amp_type='cyl', b_0=2.1, b_amp_param=8.0, pitch=11.5, rho_0=8.0, xi_0=25.0, z0=1.0):
        self.b_0 = b_0
        self.b_amp_param = b_amp_param
        self.b_amp_type = b_amp_type
        self.pitch = pitch
        self.rho_0 = rho_0
        self.xi_0 = xi_0
        self.z0 = z0
        super().__init__()

    def evaluate_grid(self, grid_x, grid_y, grid_z):
        # This will not work, still need to wrap "evaluate_function_on_grid"
        return self._evaluate_grid(grid_x, grid_y, grid_z, self.at_position)

    def at_position(self, x, y, z):
        rho, z, phi = cart2cyl([x, y, z])
        b_amp = self.b0_of_r(rho, z)
        ### Tilt angle
        xi_z = self.xi_0 * np.tanh(z / self.z0)

        # cylindrical component of the magnetic vector field
        b_rho = b_amp * np.sin(self.pitch) * np.cos(xi_z)
        b_phi = b_amp * np.cos(self.pitch) * np.cos(xi_z)
        b_z = b_amp * np.sin(xi_z)

        return [b_rho, b_z, b_phi]

    def b0_of_r(self, rho, z):

        """
        Internal function
            Intended to modulate the magnetic field amplitude by a function
            of the radial (cyl or sph) coordinate
        If constant:     B_amp = B_0     for all rho,z
        If cylindrical:  B_amp = B_0 * 1/(1 + rho/B_amp_param)
        If spherical:    B_amp = B_0 * exp(-(r - rho_0)/B_amp_param)
        B_amp = B_0    if constant
                B_0 * 1/(1 + rho/rho_0)    if cylindrical
                B_0 * exp(-(r - rho_0)/B_amp_param) if spherical
        B_amp is automatically normalized such that B_amp(at sun) = B_sun meant
        to be given by the 'B_0' param in kwargs.
        INPUT:
        ------
          - rho : cylindircal radial coordinate
          - z : height coordinates
          **kwargs : containing the information to build the wanted function
             - B_0 : an overall amplitude
             - B_amp_param : an additional parameter for the radial function
             - B_amp_type : string to specify the fonctional form. Should be
                            'cst','cyl' or 'sph'
             - rho_0 : is supposed to contained the dist. btw the Sun and the GC
        OUTPUT:
        ------
          - B_amp : the field amplitude at each location specified by rho,z
        Creation date : Jul 5 2017
        @author: V.Pelgrims
        """
        ###
        if self.b_amp_type == 'cst':
            b_amp = self.b_0
        elif self.b_amp_type == 'cyl':
            b_0 = self.b_0 * (1. + self.rho_0 / self.b_amp_param)
            b_amp = b_0 * 1. / (1. + rho / self.b_amp_param)
        elif self.b_amp_type == 'sph':
            b_amp = (self.b_0 * np.exp(-((rho ** 2 + z ** 2) ** .5 - self.rho_0) / self.b_amp_param))
        else:
            print(self.b_amp_type)
            raise ValueError('''
            Bad entry for optional argument 'B_amp_type'.
            Key must be one either: 'cst', 'sph' or 'cyl' ''')

        return b_amp
