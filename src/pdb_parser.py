#!/usr/bin/env python3

import argparse
import numpy as np


def pdb_parser(pdb_file):
    """
    This function manage pdb file isolating alpha carbons only (atome type).
    
    It retrieves:
    1. residu number
    2. x, y, z coordinates

    Args:
        pdb_file : file with pdb extension

    Returns:
        Dictionnary: {atom_number : array[coord_x, coord_y, coord_z]}
    """

    with open(pdb_file, "r") as fillin:
        coordinates_dict = {}
        for line in fillin :
            atom_type = line[12:16].strip()
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
    parser.add_argument("-p", "--pdb-file", help="file .pdb is required ", type=str)
    args = parser.parse_args()
    PDB_FILE_NAME = args.pdb_file
    if PDB_FILE_NAME[-4:] != '.pdb':
        print("Please select file with pdb extension")
    else:
        pdb_coordinates = pdb_parser(PDB_FILE_NAME)
        
