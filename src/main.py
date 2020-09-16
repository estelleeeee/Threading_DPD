import argparse

import traitement_sequences
import matrice_threading

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




liste_AA = code_1_lettre_3_lettres('seq_cible.fasta')
#print(liste_AA)

dope_energy_value = dope_parser('dope.par.txt')
pdb_file = pdb_parser('1n0a.pdb')
#print(pdb_file) #bug check
mat = matrice_distance(pdb_file)
#print(mat)
#print(len(mat[0]))
"""
for num_AA in range(len(liste_AA)):
    matrice_dope_score = mat_dope_res_pos(mat,liste_AA,liste_AA[num_AA], 1)
    print(matrice_dope_score)
"""
matrice_dope_score = mat_dope_res_pos(mat,liste_AA,"GLU", 9)
#print(matrice_dope_score)
"""
for num_AA in range(len(liste_AA)):
    for num_pos in range(1, 11):
        matrice_dope_score = mat_dope_res_pos(mat, liste_AA, liste_AA[num_AA], num_pos)
"""