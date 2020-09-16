"""This file computes the double dynamic programming

USAGE
=====
    python main.py -p <pdb_file> -d <dope_score_file> -f <fasta_file>

"""

__authors__ = ("Theo Ferreira", "Estelle Mariaux")
__date__ = "2020/09"

import argparse

import parsing
import distance
import alignment

if __name__ == "__main__":
    parser = argparse.ArgumentParser("This program aligns a sequence to a model"
                                +" structure using threading using double dynamic program."
                                +" It is based on THREADER by Jones DT 1998")
    parser.add_argument("-p", "--pdb-file",
                        help="file .pdb is required ", type=str, required = True)
    parser.add_argument("-d", "--dope-file",
                        help="dope.par.txt ", type=str, required = True)
    parser.add_argument("-f", "--fasta-file",
                        help="file .fasta is required ", type=str, required = True)
    args = parser.parse_args()

    PDB_FILE_NAME = args.pdb_file
    DOPE_FILE_NAME = args.dope_file
    FASTA_FILE_NAME = args.fasta_file

    if not PDB_FILE_NAME.endswith('.pdb'):
        print("Please select file with pdb extension")
    if not FASTA_FILE_NAME.endswith('.fasta'):
        print("Please select file with pdb extension")
    else:
        pdb_coordinates = parsing.pdb_parser(PDB_FILE_NAME)
        dope_energy_value = parsing.dope_parser(DOPE_FILE_NAME)
        seq_3_list = parsing.fasta_code(FASTA_FILE_NAME)
        mat_dist = distance.distance_matrix(pdb_coordinates)

        L_matrix = alignment.matrix_low_level(mat_dist, seq_3_list,\
                seq_3_list[0], 1, dope_energy_value)
        print(L_matrix)

        print('\nFor each amino acids, for the position 1, 4 or 7, we calculate the minimized score')
        for amino_acid in seq_3_list:
            print('-----')
            for pos in (1, 4, 7):
                alignment.matrix_low_level(mat_dist, seq_3_list, amino_acid, pos, dope_energy_value)
