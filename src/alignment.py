"""
This module creates the low level matrix for a sequence and a template.
"""

__authors__ = ("Theo Ferreira", "Estelle Mariaux")
__date__ = "2020/09"

import math
import numpy as np

import parsing

def matrix_low_level(dist_matrix, query, res_i, pos, energy_value_dict):
    '''
    This function does low level matrix for each residue in each position

    Parameters
    ----------
    dist_matrix: numpy.ndarray
        Matrix which contains distance between template coordinates
    query: dictionary
        contains query sequence
    res_i: str
        residue which was fixed
    pos: int
        position of the fixed residue

    Returns
    -------
    Low-level matrix and the optimized score for each residue in each position
    '''
    # distance between fixed residue and other positions
    dist_pos = dist_matrix[pos - 1]
    dope_mat = np.empty((len(query)+1, len(dist_pos)+1))
    dope_mat.fill(np.nan)
    dope_mat[0, :] = 0
    dope_mat[:, 0] = 0
    for row_prev in range(1, pos):
        for col_prev in range(1, pos):
            dope_mat[row_prev][col_prev] = min(
                dope_mat[row_prev-1][col_prev],
                dope_mat[row_prev][col_prev-1],
                dope_mat[row_prev-1][col_prev-1]
                    + parsing.dope_score(res_i, query[row_prev-1],\
                    dist_pos[col_prev-1], energy_value_dict))
    # residue fixed, match score = 0
    dope_mat[pos][pos] = dope_mat[pos-1][pos-1]
    for row_after in range(pos+1, len(query)+1):
        for col_after in range(pos+1, len(dist_pos)+1):
            if math.isnan(dope_mat[row_after-1][col_after])\
                    and math.isnan(dope_mat[row_after][col_after-1]):
                dope_mat[row_after][col_after] =\
                    dope_mat[row_after-1][col_after-1] +\
                    parsing.dope_score(res_i, query[row_after], \
                    dist_pos[col_after], energy_value_dict)
            elif math.isnan(dope_mat[row_after-1][col_after-1])\
                    and math.isnan(dope_mat[row_after-1][col_after]):
                dope_mat[row_after][col_after] =\
                    dope_mat[row_after][col_after-1]
            elif math.isnan(dope_mat[row_after-1][col_after-1])\
                    and math.isnan(dope_mat[row_after][col_after-1]):
                dope_mat[row_after][col_after] =\
                    dope_mat[row_after-1][col_after]
            else:
                dope_mat[row_after][col_after] = min(
                    dope_mat[row_after-1][col_after],
                    dope_mat[row_after][col_after-1],
                    dope_mat[row_after-1][col_after-1]
                        + parsing.dope_score(res_i, query[row_after-1], \
                        dist_pos[col_after-1], energy_value_dict))
    print(
        "For the fixed residue is {:s} in position {:d}, the optimized score is {:.2f}"
        .format(res_i, pos, dope_mat[-1][-1]))
    return dope_mat
