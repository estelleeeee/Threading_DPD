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