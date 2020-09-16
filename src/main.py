"""This file computes the double dynamic programming

USAGE
=====
    python main.py -p pdb_file -d ...
    

"""

__authors__ = ("Theo Ferreira", "Estelle Mariaux")
__date__ = "2020/09"

import alignment
import distance
import parsing

import sys

if __name__ == "__main__":
    pdb_coordinates = parsing.pdb_parser('1n0a.pdb')
    print(pdb_coordinates)
    #dope_energy_value = parsing.dope_parser()
    #print(dope_energy_value)
    mat_dist = distance.distance_matrix(pdb_coordinates)
    print(mat_dist)
    list_AA = parsing.fasta_code('1n09.fasta')
     
    
    matrice_dope_score = alignment.matrix_low_level(mat_dist,list_AA,"GLU", 5)
    print(matrice_dope_score)



"""
    parser = argparse.ArgumentParser(
        "Manage pdb files to isolate alpha carbons")
    parser.add_argument("-p", "--pdb-file",
                        help="file .pdb is required ", type=str)
    args = parser.parse_args()
    PDB_FILE_NAME = args.pdb_file
    if PDB_FILE_NAME[-4:] != '.pdb':
        print("Please select file with pdb extension")
    else:
        pdb_coordinates = parsing.pdb_parser(PDB_FILE_NAME)
    
    parser = argparse.ArgumentParser(
        "Manage the dope.par.txt file to isolate CA-CA pair energies")
    parser.add_argument("-d", "--dope-file", help="dope.par.txt ", type=str)
    args = parser.parse_args()
    DOPE_FILE_NAME = args.dope_file
    if DOPE_FILE_NAME[:] != 'dope.par.txt':
        print(
            + "Please select only the file 'dope.par.txt' provided by JC Gelly")
    else:
        dope_energy_value = parsing.dope_parser(DOPE_FILE_NAME)

if __name__ == "__main__":      
            
    parser = argparse.ArgumentParser("Manage pdp files and isolate alpha carbons")
    parser.add_argument("pdb_file", help = "file .pdb ", type=str)
    args = parser.parse_args()
    if args.pdb_file[-4:] != '.pdb' :         
        print("Please select file with pdb extension")
    else: 
        pdb_file = traitement_sequences.pdb_parser(args.pdb_file)
        #print(pdb_file) #bug check
        mat = matrice_threading.matrice_distance(pdb_file)
        print(mat)

#PDB
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


liste_AA = code_1_lettre_3_lettres('seq_cible.fasta')
#print(liste_AA)

dope_energy_value = dope_parser('dope.par.txt')
pdb_file = pdb_parser('1n0a.pdb')
#print(pdb_file) #bug check
mat = matrice_distance(pdb_file)
#print(mat)
#print(len(mat[0]))

for num_AA in range(len(liste_AA)):
    matrice_dope_score = mat_dope_res_pos(mat,liste_AA,liste_AA[num_AA], 1)

matrice_dope_score = mat_dope_res_pos(mat,liste_AA,"GLU", 9)
#print(matrice_dope_score)

for num_AA in range(len(liste_AA)):
    for num_pos in range(1, 11):
        matrice_dope_score = mat_dope_res_pos(mat, liste_AA, liste_AA[num_AA], num_pos)
"""