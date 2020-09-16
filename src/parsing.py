"""
This module parse pdb, dope and fasta files.
"""

__authors__ = ("Theo Ferreira", "Estelle Mariaux")
__date__ = "2020/09"


import numpy as np

DOPE_FILE = 'dope.par.txt'
AA_DICT = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
               'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
               'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
               'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

def pdb_parser(pdb_file):
    '''
    This function manage pdb file isolating alpha carbons only (atome type).

    It retrieves:
    1. residu number
    2. x, y, z coordinates

    Args:
        pdb_file : file with pdb extension

    Returns:
        Dictionnary: {atom_number : array[coord_x, coord_y, coord_z]}
    '''
    with open(pdb_file, "r") as fillin:
        coordinates_dict = {}
        for line in fillin:
            atom_type = line[12:16].strip()
            # HETATOM were excluded and only C-alpha were selected
            if line[0:6].strip() == "ATOM" and (
                    atom_type in ('CA', 'C1*')):
                res_number = int(line[22:26].strip())
                coordinates_dict[res_number] =\
                    np.array([float(line[30:38].strip()),
                              float(line[38:46].strip()),
                              float(line[46:54].strip())])
    return coordinates_dict

def dope_parser():
    '''
    This function manage dope file isolating alpha carbon - alpha carbon pair energies

    Statistical potential values :
    Contains energy value for a pair of specific atoms at a given bin distance of 0.5 Angstroms.
    First value is for a pair of atoms at a distance between 0.25 and 0.75 A.
    Bin distance reaches 15.0 A.

    It retrieves :
    1 - energy value for each residue - residue

    Returns :
        Dictionnary: {residue 1 - residue 2 : array[energy_value from 1 to end]}
    '''
    with open(DOPE_FILE, "r") as fillin:
        energy_value_dict = {}
        for line in fillin:
            words_list = line.strip().split(" ")
            res_first = words_list[0].strip()
            atom_first = words_list[1].strip()
            res_sec = words_list[2].strip()
            atom_sec = words_list[3].strip()
            # only C-alpha were selected
            if atom_first == "CA" and atom_sec == "CA":
                res_res = (res_first, res_sec)
                bin_list = [float(i) for i in words_list[4:]]
                energy_value_dict[res_res] = bin_list
        return energy_value_dict

def dope_score(res1, res2, dist_res1_res2):
    '''
    This function calculates the energy value (dope.par) for a distance between two amino acids.

    It retrieves :
    1 - energy value for a distance between two amino acids

    Args :
        res1 : residue 1
        res2 : residue 2
        dist_res1_res2 : distance between residue 1 and residue 2
        dope_file : dope score for AA

    Returns :
        Dictionnary: {residue 1 - residue 2 : interval_index}
        interval_index is the bin correponding to the distance
    '''
    dope_dict = dope_parser()
    interval_index = round(int(dist_res1_res2 * 2 - 0.5))
    return dope_dict[(res1, res2)][interval_index]


def fasta_code(fasta_file):
    '''
    This Function transforms a fasta file with AA in 1 letter code
    into a dictionary with AA, code 3 letters, in a list.

    Parameters
    ----------
    fasta_file : file.fasta

    Returns
    -------
    dict[str]=[str,str,...]
    '''
    with open(fasta_file, "r") as file:
        seq_3_list = []
        for line in file:
            if line[0] != ">":
                for amino in line[:-1]:
                    # avoid unknown residues
                    if amino != 'X':
                        seq_3_list.append(AA_DICT[amino])
    return seq_3_list
