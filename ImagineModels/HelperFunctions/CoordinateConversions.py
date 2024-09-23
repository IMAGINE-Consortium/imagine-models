import numpy as np


def cyl2cart(coordinate):
    """
    Converts cylindrical coordinates to Cartesian coordinates.
    INPUT:
    -----
     - coordinate: (3,npts)-array with cylindrical coordinates
                   [rho, phi, z] is assumed
    OUTPUT:
    -------
     - ouput : (3,npts)-array of cartesian coordinates, [x, y, z] is assumed
    @author: V.Pelgrims, adapted S.Hutschenreuter
    """
    x = coordinate[0] * np.cos(coordinate[1])
    y = coordinate[0] * np.sin(coordinate[1])
    z = coordinate[2]

    return np.array([x,y,z])


def cart2cyl(coordinate):
    """
    Converts cylindrical coordinates to Cartesian coordinates.
    INPUT:
    -----
     - coordinate: (3,npts)-array with cartesian coordinates
                   [x, y, z] is assumed
    OUTPUT:
    -------
     - ouput : (3,npts)-array of cylindrical coordinates, [rho, phi, z] is assumed
    @author: V.Pelgrims, adapted S.Hutschenreuter
    """
    rho = np.sqrt(coordinate[0]**2 + coordinate[1]**2)
    phi = np.arctan2(coordinate[1], coordinate[0]) 
    z = coordinate[2]

    return np.array([rho, phi, z])