#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module provides facilities function in order to compute various quantities
related to the discrete curvature around a given point in R^3. The discrete
curvature is defined as the angular defect around a given vertex.
"""

import numpy as np
import pandas as pd

__author__ = "Germain Salvato Vallverdu"
__copyright__ = ""
__version__ = "0.1"
__maintainer__ = "Germain Salvato Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"

units = {"angular defect": "(degrees)",
         "haddon": "(degrees)",
         "improper angle": "(degrees)",
         "distance from average plane": "(A)",
         "distances with neighbors": "(A)",
         "atom index": ""}

ATOM_COLOR_DICT = {
    "C": "#c8c8c8",
    "H": "#ffffff",
    "N": "#8f8fff",
    "S": "#ffc832",
    "O": "#f00000",
    "F": "#ffff00",
    "P": "#ffa500",
    "Cl": "#774400",
    "Br": "#774400",
    "K": "#42f4ee",
    "G": "#3f3f3f"
}

BOND_CUTOFF_DICT = {
    # distances in angstrom
    ("H", "H"): 1.0,
    ("C", "C"): 1.7,
    ("C", "O"): 1.7,
    ("C", "S"): 1.9,
    ("C", "N"): 1.7,
    ("C", "Cl"): 1.9,
    ("C", "Br"): 2.2,
    ("C", "F"): 1.7,
    ("C", "H"): 1.2,
    ("C", "P"): 1.9,
    ("C", "H"): 1.2,
    ("N", "H"): 1.2,
    ("O", "H"): 1.2,
    ("S", "H"): 1.3,
    ("P", "H"): 1.3,
    ("P", "O"): 1.8,
    ("P", "N"): 1.8,
    ("S", "O"): 1.7,
    ("F", "H"): 1.2,
}

def get_atom_color(specie):
    """
    Return the color of an atom according to the ATOM_COLOR_DICT
    """
    if specie in ATOM_COLOR_DICT:
        return ATOM_COLOR_DICT[specie]
    else:
        return ATOM_COLOR_DICT["G"]

def get_bond_cutoff(specie1, specie2, rcut=2.0):
    """
    Return the bond cut-off in angstrom corresponding to the BOND_CUTOFF_DICT
    """
    if (specie1, specie2) in BOND_CUTOFF_DICT:
        return BOND_CUTOFF_DICT[(specie1, specie2)]
    elif (specie2, specie1) in BOND_CUTOFF_DICT:
        return BOND_CUTOFF_DICT[(specie2, specie1)]
    else:
        return rcut

def read_molecule(file):
    """
    Return coords and atom names from a xyz like file. 
    The file is suposed to display the number of atom on the first line,
    followed by a title line and followed by the structure in cartesian
    coordinates. Each line contains the element as first item and the
    cartesian coordinates as 2d, 3th and 4th itmes.
    """

    natom = int(file.readline().split()[0])
    _ = file.readline()
    coords = list()
    species = list()
    for _ in range(natom):
        data = file.readline().split()[:4]
        species.append(data[0])
        coords.append(data[1:4])

    return species, np.array(coords, dtype=np.float)


def get_molecular_data(species, coords, rcut=2.0):
    """ 
    Set up the molecular data dictionnary from the species and the 
    coordinates of a molecules. Bonds are defined from a cutoff distance. If
    the bond is not available in the BOND_CUTOFF_DICT, rcut is used.
    """

    # set up json file
    model_data = {"atoms": [], "bonds": []}

    # structure part
    for iat, (specie, coord) in enumerate(zip(species, coords)):
        name = "%s%d" % (specie, iat + 1)
        model_data["atoms"].append({"name": name,
                                    "serial": iat,
                                    "element": specie,
                                    "positions": coord})

    # look for bonds
    # compute all distances and keep triangle low matrix
    coords = np.array(coords, dtype=np.float)
    distances = (coords[:, None, :] - coords[None, :, :]) ** 2
    distances = np.sqrt(np.sum(distances, axis=-1))

    natom = len(species)
    for iat in range(natom - 1):
        for jat in range(iat + 1, natom):
            rc = get_bond_cutoff(species[iat], species[jat], rcut)
            if distances[iat, jat] < rc:
                model_data["bonds"].append(
                    {"atom1_index": iat, "atom2_index": jat}
                )
    
    # distances = np.tril(np.sqrt(np.sum(distances, axis=-1)))
    # pairs = np.where((distances > 0) & (distances < rcut))
    # pairs = np.vstack(pairs).transpose()
    # model_data["bonds"] = [{"atom1_index": iat, "atom2_index": jat}
    #                        for iat, jat in pairs]

    return model_data


def compute_data(species, coords, rcut=2.0, distances=None):
    """
    Compute the data over the whole molecule. For each atom, the following data 
    are computed:
    * discrete curvature (angular defect)
    * haddon definition of curvature
    * improper angle
    * distance from the nearest plane of neighbors
    * the average distances with neighbors

    Args:
        species (list): list of elements as string
        coords (nat x 3 arry): cartesian coordinates
        rcut (float): cutoff radius to select neighbors

    """

    # compute the distance matrix
    if not distances:
        # compute all distances
        distances = (coords[:, None, :] - coords[None, :, :]) ** 2
        distances = np.sqrt(np.sum(distances, axis=-1))
    
    # dict of data to set up the dataframe
    natom = len(species)
    data = {"angular defect": list(),
            "haddon": list(),
            "improper angle": list(),
            "dist. from ave. plane": list(),
            "ave. neighb. dist.": list(),
            "neighbors": list(),
            "atom index": np.arange(1, natom + 1).tolist(),
            "species": species,
            "x": coords[:, 0],
            "y": coords[:, 1],
            "z": coords[:, 2]}
    columns = ["atom index", "species", "x", "y", "z", "angular defect", 
               "haddon", "improper angle", "dist. from ave. plane", 
               "neighbors", "ave. neighb. dist."]

    for iat in range(natom):
        sites = coords[iat, :]

        # look for neighbors
        n = 0
        ave_distance = 0
        for jat in range(natom):
            if jat == iat:
                continue

            rc = get_bond_cutoff(species[iat], species[jat], rcut)
            if distances[iat, jat] < rc:
                sites = np.row_stack((sites, coords[jat, :]))
                n += 1
                ave_distance += distances[iat, jat]

        # number of neighbors
        data["neighbors"].append(n)

        # compute data if the number of neighbors is relevant
        if n == 3:
            data["haddon"].append(get_haddon(sites)[0])
            data["improper angle"].append(
                np.abs(get_improper(sites)))
        else:
            data["haddon"].append(np.nan)
            data["improper angle"].append(np.nan)

        if n >= 3:
            data["dist. from ave. plane"].append(
                np.abs(get_pyramidalization(sites)))
            data["angular defect"].append(
                get_discrete_curvature(sites))
        else:
            data["dist. from ave. plane"].append(np.nan)
            data["angular defect"].append(np.nan)

        if n >= 1:
            data["ave. neighb. dist."].append(ave_distance / n)
        else:
            data["ave. neighb. dist."].append(np.nan)

    df = pd.DataFrame(data, columns=columns)
    df.set_index(df["atom index"], inplace=True)

    return df, distances


def get_improper(coords):
    """
    Assuming i is bonded to j, k and l, compute the improper angle between planes
    defined by (i, j, k) and (j, k, l).


                  l
                  |
                  i
                 / \
                j   k

    Args:
        coords (4 x 3 array): Coordinnates of point i, j, k and l. i must be the first

    Returns:
        improper angle (degrees)
    """

    coords = np.array(coords)
    if coords.shape != (4, 3):
        raise ValueError("The shape of the input coordinates must be (4, 3) "
                         "corresponding to 4 poinst in R^3.")

    icoords, jcoords, kcoords, lcoords = coords

    v1 = kcoords - lcoords
    v2 = jcoords - kcoords
    v3 = icoords - jcoords
    v23 = np.cross(v2, v3)  # perpendicular to j, k, l
    v12 = np.cross(v1, v2)  # perpendicular to i, j, k

    theta = np.arctan2(np.linalg.norm(v2) * np.dot(v1, v23), np.dot(v12, v23))

    return np.degrees(theta)


def get_pyramidalization(coords):
    """
    Assuming the first point is linked to the following, compute the distance
    from the average plane defined by points coords[1:, :].

    Args:
        coords (N x 3 array): Coordinnates of the in R^3. The distance is compute
            between the first one and the plane defined by the followings.

    Returns:
        distance (float)
    """

    coords = np.array(coords)
    if coords.shape[0] < 4:
        raise ValueError(
            "Input coordinates must correspond to at least 4 points in R^3.")
    if coords.shape[1] != 3:
        raise ValueError("Input coordinates must correspond to point in R^3.")

    # barycenter of the points
    G = coords[1:, :].sum(axis=0) / coords[1:, :].shape[0]

    # Cpmpute the best plane from SVD
    _, _, (vecx, vecy, u_norm) = np.linalg.svd(coords[1:, :] - G)

    return np.dot(coords[0, :] - G, u_norm)


def get_haddon(coords):
    """
    Assuming the first point is linked to the followings, compute the pyramidalization
    following the definition of Haddon et al. The angle of pyramidalization is
    compute as the angle between the vector normal to the (j, k, l) plane and the
    bonds between i and j, k and l.

    Args:
        coords (4 x 3 array): Coordinnates of the in R^3.

    Returns:
        list of float corresponding to the angles
    """

    coords = np.array(coords)
    if coords.shape != (4, 3):
        raise ValueError("The shape of the input coordinates must be (4, 3) "
                         "corresponding to 4 poinst in R^3.")

    # barycenter of the points
    G = coords[1:, :].sum(axis=0) / coords[1:, :].shape[0]

    # Cpmpute the best plane from SVD
    _, _, (vecx, vecy, u_norm) = np.linalg.svd(coords[1:, :] - G)

    # change the
    GI = coords[0, :] - G
    if np.dot(GI, u_norm) < 0:
        u_norm = - u_norm

    angles = list()
    for coord in coords[1:]:
        v = coord - coords[0, :]
        v /= np.linalg.norm(v)
        cos = np.dot(u_norm, v)
        angles.append(np.degrees(np.arccos(cos)))

    haddon = np.array(angles).mean() - 90.0

    return haddon, angles


def get_discrete_curvature(coords):
    """
    Compute the discrete curvature (angular defect) around a vertex of a set of
    points in R^3. The function works for any number of points greater than 4:
    one vertex and 3 points surrounding the vertex defining the edges. The needed
    data are the list of the coordinates of the points, the first one being the
    vertex's one.
    The function first looks for the best fit plane of the points surrounding the
    vertex and sorts that points in order to compute the angles between the edges
    connected to the vertex and the the angular defect.

    Args:
        coords (ndarray N x 3): Coordinates of the points. The first one is the vertex.

    Returns
        curvature (float)
    """
    coords = np.array(coords)

    if coords.shape[1] != 3:
        raise ValueError("3D coordinates are needed.")

    npts = coords.shape[0]
    if npts < 4:
        raise ValueError("Not enough points to compute curvature")

    # barycenter of the points
    G = coords[1:, :].sum(axis=0) / coords[1:, :].shape[0]

    # Cpmpute the best plane from SVD
    _, _, (vecx, vecy, u_norm) = np.linalg.svd(coords[1:, :] - G)

    # compute all angles with vectorx
    angles = list()
    for points in coords[1:, :]:
        u = points - G
        u /= np.linalg.norm(u)

        cos = np.dot(vecx, u)
        if np.dot(vecy, u) > 0:
            angles.append(np.degrees(np.arccos(cos)))
        else:
            angles.append(360 - np.degrees(np.arccos(cos)))

    # sort points according to angles
    idx = np.arange(1, npts)
    idx = idx[np.argsort(angles)]
    idx = np.append(idx, idx[0])

    # compute curvature
    curvature = 360
    for i, j in np.column_stack([idx[:-1], idx[1:]]):
        u = coords[i, :] - coords[0, :]
        u /= np.linalg.norm(u)

        v = coords[j, :] - coords[0, :]
        v /= np.linalg.norm(v)

        cos = np.dot(u, v)
        curvature -= np.degrees(np.arccos(cos))

    return curvature
