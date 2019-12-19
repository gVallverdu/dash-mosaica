# coding: utf-8


"""
This is the core module that provides the functions who computes the
geometrical descriptors associated to the curvature.

Notations
---------


A
*(A)
point B belongs to *(A)
I : barycenter of points in *(A)

"""

import warnings
import numpy as np

__author__ = "Germain Salvato-Vallverdu"
__copyright__ = ""
__version__ = "0.1"
__maintainer__ = ""
__email__ = "germain.vallverdu@univ-pau.fr"
__status__ = "Development"
__date__ = "June, 2018"

__all__ = ["get_improper", "get_pyr_angle", "get_pyr_distance",
           "get_angular_defect", "get_spherical_curvature",
           "get_hybridization_coeff", "get_hybridization_numbers",
           "hybridization"]


def check_array(a, shape, dtype=np.float):
    """Check array a and return a numpy array."

    Args:
        a: iterable to convert in numpy array
        shape (tuple of int): if one value is negative, it is considered as a
            min condition.
        dtype: a numpy type, default is `np.float`
    """
    # convert in numpy array
    try:
        a = np.asarray(a, dtype=dtype)
    except ValueError:
        raise ValueError("Cannot convert {} in a numpy array of float.".format(a))

    # check shape
    if np.all(np.array(shape) > 0):
        if a.shape != shape:
            raise ValueError("The shape of coordinates is {}. The expected shape is {}.".format(a.shape, shape))
    else:
        for idim, (dim, dim_a) in enumerate(zip(shape, a.shape), 1):
            if dim < 0 and dim_a < -dim:
                raise ValueError("The shape of coordinates is {}. The length in the dimension {} must be greater or equal to {}.".format(a.shape, idim, -dim))
            if 0 < dim != dim_a:
                raise ValueError("The shape of coordinates is {}. The length in dimension {} must be equal to {}.".format(a.shape, idim, dim))

    return a


def center_of_mass(coords, masses=None):
    """Compute the center of mass of the points at coordinates `coords` with
    massess `masses`.

    Args:
        coords (np.ndarray): (N, 3) matrix of the points in R^3
        masses (np.ndarray): vector of length N with the masses

    Returns:
        The center of mass, a vector in R^3
    """
    if masses is None:
        masses = np.ones(coords.shape[0])
    return np.sum(coords * masses[:, np.newaxis], axis=0) / masses.sum()


def circum_center(coords):
    """Compute the coordinates of the center of the circumscribed circle from
    three points A, B and C in R^3.
    
    Args:
        coords (ndarray): (3x3) cartesian coordinates of points A, B and C. 
            One point per line.

    Returns
        The coordinates of the center of the cicumscribed circle
    """
    a, b, c = coords

    ABvAC = np.cross(b - a, c - a)
    M = np.array([b - a, c - a, ABvAC])
    B = np.array([np.dot(b - a, (b + a) / 2),
                  np.dot(c - a, (c + a) / 2),
                  np.dot(ABvAC, a)])

    return np.dot(np.linalg.inv(M), B)


def get_plane(coords, masses=None):
    """Given a set of N points in 3D space, compute an orthonormal basis of 
    vectors, the first two belonging to the plane and the third one being normal 
    to the plane. In the particular case where N equal 3, there is an exact
    definition of the plane as the three points define an unique plan.
    
    If N = 3, use a gram-schmidt orthonormalization to compute the vectors. If
    N > 3, the orthonormal basis is obtained from SVD.

    Args:
        coords (np.ndarray): (N, 3) matrix of the points in R^3
        masses (np.ndarray): vector of length N with the masses
    
    Returns:
        The function returns the orthonormal basis (vecx, vecy, n_a), vector
        n_a being normal to the plane.
    """
    if masses is None:
        masses = np.ones(coords.shape[0])
    com = center_of_mass(coords, masses)

    if coords.shape == (3, 3):
        # the plane is exactly defined from 3 points
        vecx = coords[1] - coords[0]
        vecx /= np.linalg.norm(vecx)

        # vecy, orthonormal with vecx
        vecy = coords[2] - coords[0]
        vecy -= np.dot(vecy, vecx) * vecx
        vecy /= np.linalg.norm(vecy)

        # normal vector
        n_a = np.cross(vecx, vecy)

    else:
        # get the best fitting plane from SVD.
        _, _, (vecx, vecy, n_a) = np.linalg.svd(coords - com)

    return vecx, vecy, n_a


def regularize(a, star_a):
    """Regularize the coordinates of the points in *(A) such as all distances 
    between points in *(A) and A are equal.
    
    Args:
        a (np.ndarray): Cartesian coordinates of atom a
        star_a (nd.array): Cartesian coordinates of atoms in *(A)

    Returns
        new coordinates of star_a
    """

    u = star_a - a
    norm = np.linalg.norm(u, axis=1)
    u /= norm[:, np.newaxis]

    return a + u


def get_central_atom(coords):
    """Try to guess the central atom. The central atom is the one which is
    bonded to all other atoms. Here we assume that the central atom is the
    closest to the center of mass.

    WARNINGS: This condition is not completely sufficient and can lead to wrong
    results in certains particular case.

    Args:
        coords (np.ndarray): cartesian coordinates in R^3, shape (N x 3)

    Returns:
        The index of the central atom.
    """
    warnings.warn(("This function *tries* to return the central point from"
                   " a given set of point but it may fail in some special "
                   "case or specific geometries. Be careful!"),
                  category=UserWarning)

    com = center_of_mass(coords)
    distances = np.sum((coords - com[np.newaxis, :])**2, axis=1)
    return np.argmin(distances)
    # distances = (coords[np.newaxis, :, :] - coords[:, np.newaxis, :])**2
    # distances = np.sqrt(np.sum(distances, axis=2))
    #
    # return np.argmin(distances.sum(axis=1))


def get_improper(a, star_a):
    """Compute the improper angle between planes defined by atoms (i, j, k) and
    (j, k, l). Atom A, is atom i and atoms i, j and k belong to *(A).

                  l
                  |
                  i
                 / \
               j     k

    Args:
        a (np.ndarray): Cartesian coordinates of atom a
        star_a (nd.array): (3x3) matrix of cartesian coordinates of atoms i, j, k

    Returns:
        improper angle (degrees)
    """
    icoords = check_array(a, shape=(3,), dtype=np.float)
    star_a = check_array(star_a, shape=(3, 3), dtype=np.float)

    star_a = regularize(a, star_a)

    # if reorder:
    #     # get central atom and put central atoms as first atom
    #     iat = get_central_atom(coords)
    #     coords[[0, iat]] = coords[[iat, 0]]

    jcoords, kcoords, lcoords = star_a

    # compute vectors
    vij = jcoords - icoords
    vjk = kcoords - jcoords
    vlk = kcoords - lcoords
    m = np.cross(vij, vjk)  # perpendicular to ijk
    n = np.cross(vlk, vjk)  # perpendicular to jkl

    # print("vij ", vij)
    # print("vjk ", vjk)
    # print("vlk ", vlk)
    # print("m   ", m)
    # print("n   ", n)

    # compute the angle
    theta = np.arctan2(np.dot(vij, n) * np.linalg.norm(vjk), np.dot(m, n))
    # theta2 = np.arccos(np.dot(m, n) / np.linalg.norm(m) / np.linalg.norm(n))
    # print(np.degrees(theta), np.degrees(theta2))

    return np.degrees(theta)


def get_pyr_distance(a, star_a):
    """Compute the distance between atom A and the plane define by *(A) or
    the best fitting plane of *(A).

    Args:
        a (np.ndarray): Cartesian coordinates of atom a
        star_a (nd.array): (N x 3) cartesian coordinates of atoms in *(A)

    Returns:
        distance (float) in the same unit as the input coordinates.
    """
    a = check_array(a, shape=(3,), dtype=np.float)
    star_a = check_array(star_a, shape=(-3, 3), dtype=np.float)

    _, _, n_a = get_plane(star_a)

    return np.abs(np.dot(a - center_of_mass(star_a), n_a))


def get_pyr_angle(a, star_a):
    """Compute the pyramidalization angle considering a pyramidal structure
    defined from a set of points in R^3. Point A is the vertex of
    the pyramid, other points belong to *(A).

    Args:
        a (np.ndarray): Cartesian coordinates of atom a
        star_a (nd.array): (N x 3) cartesian coordinates of atoms in *(A)

    Returns:
        The pyramidalization angle in degrees
    """
    # check input coords
    a = check_array(a, dtype=np.float, shape=(3,))
    star_a = check_array(star_a, dtype=np.float, shape=(-3, 3))

    # regularize coords
    star_a = regularize(a, star_a)

    # get the normal vector to the plane defined from *(A)
    _, _, n_a = get_plane(star_a)

    # change the direction of n_a to be the same as IA.
    IA = a - center_of_mass(star_a)
    n_a = -n_a if np.dot(IA, n_a) < 0 else n_a

    v = star_a[0] - a
    v /= np.linalg.norm(v)
    theta = np.degrees(np.arccos(np.dot(v, n_a)))

    return theta - 90.


def get_angular_defect(a, star_a):
    """
    Compute the angular defect as a measure of the discrete curvature around a
    vertex, point A in R^3, and points B in R^3 belonging to *(A) and bonded to
    A.
    The function first looks for the best fit plane of points in *(A)
    and sorts that points in order to compute the angles between the edges
    connected to the vertex (A).

    Args:
        a (np.ndarray): Cartesian coordinates of atom a
        star_a (nd.array): (N x 3) cartesian coordinates of atoms in *(A)

    Returns
        curvature (float)
    """
    # check and regularize coords
    a = check_array(a, shape=(3,), dtype=np.float)
    star_a = check_array(star_a, shape=(-3, 3), dtype=np.float)
    star_a = regularize(a, star_a)

    # get P the plane of *(A)
    vecx, vecy, _ = get_plane(star_a)

    # compute all angles with vecx in order to sort atoms of *(A)
    u = star_a - center_of_mass(star_a)
    norm = np.linalg.norm(u, axis=1)
    u /= norm[:, np.newaxis]
    cos = np.dot(u, vecx)
    angles = np.where(np.dot(u, vecy) > 0,
                      np.degrees(np.arccos(cos)),
                      360 - np.degrees(np.arccos(cos)))

    # sort points according to angles
    idx = np.arange(angles.size)
    idx = idx[np.argsort(angles)]
    idx = np.append(idx, idx[0])

    # compute curvature
    ang_defect = 360
    for i, j in np.column_stack([idx[:-1], idx[1:]]):
        u = star_a[i, :] - a
        u /= np.linalg.norm(u)

        v = star_a[j, :] - a
        v /= np.linalg.norm(v)

        cos = np.dot(u, v)
        ang_defect-= np.degrees(np.arccos(cos))

    return ang_defect


def get_spherical_curvature(a, star_a):
    """Compute the spherical curvature associated to the osculating sphere at
    a give point A of a molecule bonded to atoms B belonging to *(A).

    Args:
        a (np.ndarray): Cartesian coordinates of atom a
        star_a (nd.array): (3 x 3) cartesian coordinates of atoms in *(A)

    Returns
        curvature (float)
    """
    # check and regularize coords
    a = check_array(a, shape=(3,), dtype=np.float)
    star_a = check_array(star_a, shape=(3, 3), dtype=np.float)

    # plane *(A)
    point_O = circum_center(star_a)
    _, _, n_a = get_plane(star_a)

    # needed length
    l = np.linalg.norm(star_a[0] - point_O)
    z_A = np.dot(a - point_O, n_a)
    L = np.sqrt(np.linalg.norm(a - point_O)**2 - z_A**2)

    # spherical curvature
    kappa = 1 / np.sqrt(l**2 + (L**2 + z_A**2 - l**2)**2 / (4 * z_A**2))

    return kappa


def get_hybridization_coeff(pyrA=None, a=None, star_a=None, radians=False):
    """Compute the hybridization coefficient of the h_pi hybrid orbital that
    is directed along the POAV vector, normal to *(A).

    Args:
        pyrA (float): the pyramidalization angle, check radians for the unit.
        a (np.ndarray): Cartesian coordinates of atom a
        star_a (nd.array): (3 x 3) cartesian coordinates of atoms in *(A)
        radians (bool): if True, pyrA is given in radians.
    
    Returns:
        c_pi and lambda_pi the coefficient of the s and p_z atomic orbitals
        respectively.
    """
    if pyrA is not None:
        # compute hybridization coeff from pyrA angle
        r_pyrA = pyrA if radians else np.radians(pyrA)

        # check domain definition of lambda_pi
        if np.tan(r_pyrA) > 1 / np.sqrt(2):
            raise ValueError("pyrA = {} degrees lambda_pi is not define.".format(np.degrees(r_pyrA)))

        c_pi = np.sqrt(2) * np.tan(r_pyrA)
        lambda_pi = np.sqrt(1 - 2 * np.tan(r_pyrA) ** 2)
        
        return c_pi, lambda_pi

    elif a is not None and star_a is not None:
        # compute first pyramidalization angle and the hybridization coefficients
        pyrA = get_pyr_angle(a, star_a)

        return get_hybridization_coeff(pyrA)
    
    else:
        raise ValueError("You have to provide either pyrA or both a and star_a.")
    

def get_hybridization_numbers(pyrA=None, a=None, star_a=None, radians=False):
    """Compute the hybridization numbers m and n corresponding to atom A. The
    number m is linked to the hybridization coefficients and n = 3m + 2.

    Args:
        pyrA (float): the pyramidalization angle, check radians for the unit.
        a (np.ndarray): Cartesian coordinates of atom a
        star_a (nd.array): (3 x 3) cartesian coordinates of atoms in *(A)
        radians (bool): if True, pyrA is given in radians.
    
    Returns:
        m and n
    """
    if pyrA is not None:
        # compute hybridization coeff from pyrA angle
        c_pi, lambda_pi = get_hybridization_coeff(pyrA, radians=radians)

        m = (c_pi / lambda_pi) ** 2
        n = 3 * m + 2
        return m, n

    elif a is not None and star_a is not None:
        # compute first pyramidalization angle and then hybridization numbers
        pyrA = get_pyr_angle(a, star_a)

        return get_hybridization_numbers(pyrA)
    
    else:
        raise ValueError("You have to provide either pyrA or both a and star_a.")


def hybridization(pyrA=None, a=None, star_a=None, radians=False):
    """Compute the hybridization of one atom.

    Args:
        pyrA (float): the pyramidalization angle, check radians for the unit.
        a (np.ndarray): Cartesian coordinates of atom a
        star_a (nd.array): (3 x 3) cartesian coordinates of atoms in *(A)
        radians (bool): if True, pyrA is given in radians.

        s^{1-C_pi^2} p^{n-3m+C_pi^2} 
        s p^{(n - 3 m + C_pi^2) / (1 - C_pi^2)}
    
    Returns:
        h the coefficient of p in (s p^h)
    """
    if pyrA is not None:
        # compute hybridization coeff from pyrA angle
        c_pi, _ = get_hybridization_coeff(pyrA, radians=radians)
        m, n = get_hybridization_numbers(pyrA, radians=radians)

        p_coeff = (n - 3 * m + c_pi ** 2) / (1 - c_pi ** 2)
        return p_coeff

    elif a is not None and star_a is not None:
        # compute first pyramidalization angle and then hybridization numbers
        pyrA = get_pyr_angle(a, star_a)

        return get_hybridization_numbers(pyrA)
    
    else:
        raise ValueError("You have to provide either pyrA or both a and star_a.")
