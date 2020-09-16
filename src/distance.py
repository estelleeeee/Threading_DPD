"""
This module creates the distance matrix.
"""

__authors__ = ("Theo Ferreira", "Estelle Mariaux")
__date__ = "2020/09"


import math
import numpy as np


def distance_matrix(amino_dict):
    '''
    This function creates a distance matrix between all amino acids in the dictionary

    Parameters
    ----------
    AA_dict : dictionary
        contains position and coordinates of amino acids

    Returns
    -------
    np.array(dist_matrix)
        Matrix which contains distance between amino acids
    '''
    dist_matrix = []
    for num1 in amino_dict:
        dist_line = []
        for num2 in amino_dict:
            dist_line.append(distance_array(amino_dict[num1], amino_dict[num2]))
        dist_matrix.append(dist_line)
    return np.array(dist_matrix)

def distance_array(array1, array2):
    '''
    This function calculates the distance beween two arrays

    Parameters
    ----------
    array1 : numpy.ndarray
        coordinates of the first position
    array2 : numpy.ndarray
        coordinates of the second position

    Returns
    -------
    np.array(dist)
        distance between the two arrays
    '''
    x_coord = array1[0] - array2[0]
    y_coord = array1[1] - array2[1]
    z_coord = array1[2] - array2[2]
    dist = math.sqrt(x_coord**2 + y_coord**2 + z_coord**2)
    return np.array(dist)
