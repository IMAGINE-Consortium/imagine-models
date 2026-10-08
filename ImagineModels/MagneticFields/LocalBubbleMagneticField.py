import healpy as hp
import numpy as np

from ImagineModels import RegularVectorField


class LBMagneticField(RegularVectorField):
    r"""
        ===  LBMagneticField : GMF in the shell of the LocalBubble ===

        Implementation of the analytical model for the GMF in the shell of the
            Local Bubble as presented in Alves et al. 2018
                and updated in Pelgrims et al. 2020 for any close shape of the
                inner surface of the Local Bubble.
            Refs:   - https://ui.adsabs.harvard.edu/abs/2018A%26A...611L...5A/abstract
                    - https://ui.adsabs.harvard.edu/abs/2020A%26A...636A..17P/abstract

            LBMagneticField(SkyCoordinates,**kwargs{numerous parameters})

                " B_today \propto n x (B0 x er)
                "           = B0 (n.er) - er (n.B0)
                "   where   B_today = vector field now
                "           B0 is vector field before explosion (in Cart)
                "           n is normal vector to the surface of today's
                "               shell (in Cart)
                "           er is the radial vector from the explosion center
                "               to the surface point (in Cart)


            INPUT:
            ------
              [model parameter]
               - dx,dy,dz: offset of the explosion center from the Sun [kpc]
               - theta_B0,phi_B0: co-latitude and longitude describing the 3D
                                    orientation of the Bfield BEFORE explosion.
                                    (assumed to be uniform) [radians]
               - LB_distance: path to the fits file with healpix maps of the
                                heliocentric distances of the inner of the LB
                                as measured from the Sun and as determined
                                in Pelgrims+2020. [! the maps are in pc]
                             or healpix [ring format] map of wanted shape [kpc]
               - surfaceLabel: if 'path' is given, the label is used to select
                                a given model of the shape of the inner shell
                                of the LB.
                            surfaceLabel can be lmax2,lmax4,lmax6,lmax8,lmax10
                                to point to the corresponding shell model up to
                                the maximum multipole sph. harm expansion as
                                defined in Pelgrims+2020.
                    Note: the Pelgrims+2020 models for the shape of the shell
                    -----   can be downloaded from:
                                https://doi.org/10.7910/DVN/RHPVNC

            Created on Apr 11 2023
            @author: V.Pelgrims
            """

    def __init__(self, LB_distance,\
                        theta_B0=1.278, phi_B0=1.278,\
                            dx=0.058,dy=0.079,dz=-0.086,\
                             surfaceLabel='lmax6'):
        super().__init__()
        self.theta_B0 = theta_B0    # [rad]
        self.phi_B0 = phi_B0        # [rad]
        self.dx = dx                # [kpc]
        self.dy = dy                # [kpc]
        self.dz = dz                # [kpc]

        if type(LB_distance) is str:
            if surfaceLabel == 'lmax2':
                field = 1
            elif surfaceLabel == 'lmax4':
                field = 2
            elif surfaceLabel == 'lmax6':
                field = 3
            elif surfaceLabel == 'lmax8':
                field = 4
            elif surfaceLabel == 'lmax10':
                field = 5
            else:
                raise ValueError('un-recognized label for the sought LB shell map')
            #
            r_edge = hp.read_map(LB_distance,field=field) / 1000   # to be in [kpc]
            #
        elif type(LB_distance) is np.ndarray:
            r_edge = LB_distance
        else:
            raise ValueError('Bad intry for LB_surface variable.',\
                                'We expect a path to Pelgrims+2020 map',\
                                    'or an healpix map with distance to the',\
                                        'shell of your bubble')

        self.nside = hp.get_nside(r_edge)
        self.npix = hp.nside2npix(self.nside)

        # call the function that updates the surface and all related quantities
        self.UpdateShellModel(r_edge)

        #
    #

    def at_LonLat(self,SkyCoordinates=None):
        '''
            Evaluates the model toward Galactic Longitudes and Latitudes
                given in SkyCoordinates (2,n)
                Longitudes and Latitudes are expected to be in degrees.
                If not specified, full-sky is assumed.

            OUTPUT:
            -------
            [Bx,By,Bz]: the component of the GMF vector field in the shell of the LB
                            in CARTESIAN coordinates centered on the Galactic Center

        '''
        # 1. find the indices of the corresponding healpix pixels
        if SkyCoordinates is not None:
            pix_ids = hp.ang2pix(self.nside,SkyCoordinates[0],SkyCoordinates[1],\
                                    lonlat=True)
        else:
            # if sky coord. not specified, full-sky is assumed.
            pix_ids = np.arange(self.npix)
        #

        # 2. evaluates the model for those pixels

        # surface xyz in explosion center
        xell = self.Surf_Cart[0] - self.dx
        yell = self.Surf_Cart[1] - self.dy
        zell = self.Surf_Cart[2] - self.dz
        # sph coord in explosion centered ref. frame
        tell = np.mod(np.arctan2((xell**2 + yell**2)**.5,zell),np.pi)
        pell = np.mod(np.arctan2(yell,xell),2*np.pi)

        # computation of the radial unit vector from the explosion center
        erell_x = np.sin(tell)*np.cos(pell)
        erell_y = np.sin(tell)*np.sin(pell)
        erell_z = np.cos(tell)

        # the normal to the surface is given as an entry
        nell_x = self.normal[0]
        nell_y = self.normal[1]
        nell_z = self.normal[2]

        # initial B field in Cart. coord
        B0x = np.sin(self.theta_B0)*np.cos(self.phi_B0)
        B0y = np.sin(self.theta_B0)*np.sin(self.phi_B0)
        B0z = np.cos(self.theta_B0)

        # computation of necessary dot products:
        # nell . B0
        nell_B0 = nell_x*B0x + nell_y*B0y + nell_z*B0z
        # nell . erell
        nell_erell = nell_x*erell_x + nell_y*erell_y + nell_z*erell_z

        # Cut the surface and the normal and derived vectors: n and erell
        nell_x = nell_x[pix_ids]
        nell_y = nell_y[pix_ids]
        nell_z = nell_z[pix_ids]
        erell_x = erell_x[pix_ids]
        erell_y = erell_y[pix_ids]
        erell_z = erell_z[pix_ids]
        nell_erell = nell_erell[pix_ids]
        nell_B0 = nell_B0[pix_ids]

        # B_today = B0 (n.er) - er (n.B0)

        Bx = B0x * nell_erell - erell_x * nell_B0
        By = B0y * nell_erell - erell_y * nell_B0
        Bz = B0z * nell_erell - erell_z * nell_B0
        # some normalization
        Bn = (Bx**2+By**2+Bz**2)**.5
        Bx /= Bn
        By /= Bn
        Bz /= Bn
        # and that's it.

        return [Bx,By,Bz]


    def UpdateShellModel(self,distance2shell):
        '''
            Update the shell model and related quantities:
                - Cart. coordinates of the shell
                - normal vectors
        '''

        npix = self.npix
        nside = self.nside

        # a function to get the normal vectors
        def get_normal(modeled_edge):
            neigh = hp.get_all_neighbours(nside,np.arange(npix))[np.arange(0,8,2),:]
            xyz_ = np.asarray(hp.pix2vec(nside,np.arange(npix))*modeled_edge)
            xyz_neigh = xyz_[:,neigh]
            v_swne = xyz_neigh[:,2,:] - xyz_neigh[:,0,:]
            v_senw = xyz_neigh[:,1,:] - xyz_neigh[:,3,:]
            n = np.cross(v_swne.T,v_senw.T).T
            n /= np.sqrt(np.sum(n**2,axis=0))
            return n

        self.surface = distance2shell       # [kpc]

        # 3D Cartesian coordinates at the edge of the LB
        Surf_Cart = np.asarray(hp.pix2vec(nside,np.arange(npix))*distance2shell)

        # to draw the normal vectors outward from the observer, we need the
        # unit radial vectors
        erS,_,__ = u_sph(Surf_Cart)

        # normal vectors
        Surf_normal = get_normal(self.surface)

        # drawn outwards from the observer.
        Surf_normal *= np.sign(np.sum(erS*Surf_normal,axis=0))

        self.normal = Surf_normal
        self.Surf_Cart = Surf_Cart
        #
    #

    def position_at_LonLat(self,SkyCoordinates=None):
        '''
            Find the Cartesian coordinates where the sky and the model are
            defined for a specific set of lines of sight (sky coordinates).
        '''
        # 1. find the indices of the corresponding healpix pixels
        if SkyCoordinates is not None:
            pix_ids = hp.ang2pix(self.nside,SkyCoordinates[0],SkyCoordinates[1],\
                                    lonlat=True)
        else:
            # if sky coord. not specified, full-sky is assumed.
            pix_ids = np.arange(self.npix)
        #
        return self.Surf_Cart[:,pix_ids]

    def at_position(self, x, y, z):
        mess = '''You shoud not query this model from Cartesian coordinates.
                Instead, use the at_LonLat() function and specify the sky
                directions you are interested in by the model.
                You may also use the position_at_LonLat() to know what are
                    the Cartesian coordinates where the model is defined.'''

        raise ValueError(mess)
        return


def vec_Cart2Sph(vectorField,Position):
    '''
        Convert the components of a vector field given at Cart position
        to component of the vector field in heliocentric spherical coord.
        system
    '''

    # unit basis vector of sph coord. at position
    erS,etS,epS = u_sph(Position)

    vx = vectorField[0]
    vy = vectorField[1]
    vz = vectorField[2]

    vr = vx * erS[0] + vy * erS[1] + vz * erS[2]
    vt = vx * etS[0] + vy * etS[1] + vz * etS[2]
    vp = vx * epS[0] + vy * epS[1] + vz * epS[2]

    #

    return [vr,vt,vp]


def u_sph(coord):
    """
    For input cartesian coordinates, returns the 3-basis vectors of the
    spherical coordinate system.

    INPUT:
    ------
    - coord : an (3, N)-array of the coordinates at which the vectors
            have to be computed.

    OUTPUT:
    -------
    - u_r, u_theta, u_phi : three (3, N)-array containing the basis vectors

    Created on Jun 22 2016
    @author: V.Pelgrims

    # # #       STOLLEN from GalaxyBasics in gpempy     # # #

    """

    theta = np.arctan2((coord[0]**2 + coord[1]**2)**.5,coord[2])
    phi = np.arctan2(coord[1],coord[0])

    #compute the stuff
    Cphi = np.cos(phi)
    Sphi = np.sin(phi)
    Ctheta = np.cos(theta)
    Stheta = np.sin(theta)

    u_r = np.asarray([Cphi*Stheta, Sphi*Stheta, Ctheta])

    u_theta = np.asarray([Cphi*Ctheta, Sphi*Ctheta, -Stheta])

    u_phi = np.asarray([-Sphi, Cphi, np.zeros(Cphi.size)])

    return u_r, u_theta, u_phi
#



# That's it folk!
