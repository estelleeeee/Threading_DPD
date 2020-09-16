#!/usr/bin/env python3

import argparse
import numpy as np


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
                    atom_type == "CA" or atom_type == "C1*"):
                res_number = int(line[22:26].strip())
                coordinates[res_number] =\
                    np.array([float(line[30:38].strip()),
                              float(line[38:46].strip()),
                              float(line[46:54].strip())])
    return(coordinates_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Manage pdb files to isolate alpha carbons")
    parser.add_argument("-p", "--pdb-file",
                        help="file .pdb is required ", type=str)
    args = parser.parse_args()
    PDB_FILE_NAME = args.pdb_file
    if PDB_FILE_NAME[-4:] != '.pdb':
        print("Please select file with pdb extension")
    else:
        pdb_coordinates = pdb_parser(PDB_FILE_NAME)


def dope_parser(dope_file):
    '''
    This function manage dope file isolating alpha carbon - alpha carbon pair energies

    Statistical potential values : 
    Contains energy value for a pair of specific atoms at a given bin distance of 0.5 Angstroms.
    First value is for a pair of atoms at a distance between 0.25 and 0.75 A. 
    Bin distance reaches 15.0 A.

    It retrieves : 
    1 - energy value for each residue - residue

    Args : 
        dope_file : Statistical potential file `dope.par` given by Jean-Christophe Gelly
        extension .txt

    Returns : 
        Dictionnary: {residue 1 - residue 2 : array[energy_value from 1 to end]}
    '''


with open(dope_file, "r") as fillin:
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
    return(energy_value_dict)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Manage the dope.par.txt file to isolate CA-CA pair energies")
    parser.add_argument("-d", "--dope-file", help="dope.par.txt ", type=str)
    args = parser.parse_args()
    DOPE_FILE_NAME = args.dope_file
    if DOPE_FILE_NAME[:] != 'dope.par.txt':
        print(
            + "Please select only the file 'dope.par.txt' provided by JC Gelly")
    else:
        dope_energy_value = dope_parser(DOPE_FILE_NAME)


def dope_score(res1, res2, dist_res1_res2):
    '''
    This function calculates the energy value (dope.par) for a distance between two amino acids.

    It retrieves : 
    1 - energy value for a distance between two amino acids

    Args : 
        res1 : residue 1
        res2 : residue 2
        dist_res1_res2 : distance between residue 1 and residue 2 

    Returns : 
        Dictionnary: {residue 1 - residue 2 : interval_index}
        interval_index is the bin correponding to the distance
    '''

    dope_dict = dope_parser(dope_file)
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

    AA_DICT = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
               'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
               'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
               'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}
    AA_3_dict = {}
    with open(fasta_file, "r") as file:
        seq_3_list = []
        for line in file:
            if line[0] == ">":
                # dictionnary key is the name of the sequence
                seq_name = line[1:5]
            else:
                for AA in line[:-1]:
                    if AA != 'X':  # HETATOM were not selected
                        seq_3_list.append(AA_DICT[AA])
    AA_3_dict[seq_name] = seq_3_list
    return AA_3_dict
