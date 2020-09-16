"""This file computes the double dynamic programming

USAGE
=====
    python main.py -p pdb_file -d ...
    

"""

__authors__ = ("Theo Ferreira", "Estelle Mariaux")
__date__ = "2020/09"

from parsing import *
from distance import *
from alignment import *


if __name__ == "__main__":                  
    parser = argparse.ArgumentParser("This program aligns a sequence to a model \
                                structure using threading using double dynamic program. \
                                It is based on THREADER by Jones DT 1998")
    parser.add_argument("-p", "--pdb-file",
                        help="file .pdb is required ", type=str)
    parser.add_argument("-d", "--dope-file", 
                        help="dope.par.txt ", type=str)
    parser.add_argument("-f", "--fasta-file",
                        help="file .fasta is required ", type=str)
    args = parser.parse_args()

    PDB_FILE_NAME = args.pdb_file
    DOPE_FILE_NAME = args.dope_file
    FASTA_FILE_NAME = args.fasta_file

    if PDB_FILE_NAME[-4:] != '.pdb':
        print("Please select file with pdb extension")
    if DOPE_FILE_NAME[:] != 'dope.par.txt':
        print("Please select the file 'dope.par.txt' provided by JC Gelly")
    if FASTA_FILE_NAME[-6] != '.fasta' :
        print("Please select file with pdb extension")
    else:
        pdb_coordinates = pdb_parser(PDB_FILE_NAME)
        dope_energy_value = dope_parser(DOPE_FILE_NAME)
        seq_3_list = fasta_code(FASTA_FILE_NAME)
        mat_dist = distance_matrix(pdb_coordinates)

        print(matrix_low_level(distance_matrix, seq_3_list, seq_3_list[0], 1))
        # dope_mat_matrix = []


        # for res in pdb_coordinates.keys():
        #     for j in enumerate(seq_3_list):
        #         dope_mat = matrix_low_level(distance_matrix, seq_3_list, j, j+1)





