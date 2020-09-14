#!/usr/bin/env python3

'''
Question pour lundi : 
- les Class GitHub Adam
- return
- if __name__
- def __init__
- print() à chaque ligne sans exagération ++
'''



import argparse
import numpy as np
def pdb_parser(pdb_file) :
    """
    This function manage pdb file isolating alpha carbons only.
    It retrieves:
    1. atom number
    2. x, y, z coordinates

    Args:
        pdb_file : file with pdb extension

    Returns:
        Dictionnary: {atom_number, atom_name, coord_x, coord_y, coord_z}
    """
    #print("je suis dans la fonction")

    with open(pdb_file, "r") as fillin:
        atom_number = {}
        coordinates = {}
        #print(coordinates) #bug check = ok, empty dict

        '''
        #line_start = {}
        #atom_type = {}
        #res_type = {}
        #chains = {}
        #res_number = {}
        #occupency = {}
        #temperature_factors = {}
        #element_symbol = {}
        #atom_charge = {}
        '''
        
        lines = fillin.readlines()
        for line in lines:
            #print(line[0:6].strip()) #bug check = ok
            atom_type = line[12:16].strip()
            #print(atom_type) #bug check = ok
            if line[0:6].strip() == "ATOM" and (atom_type == "CA" or atom_type == "C1*"):
                #print("on est dans la condition") #bug check = ok
                atom_number = int(line[6:11].strip())
                #print(atom_number) #bug check = ok
                #line_start[atom_number] = line[0:6].strip()
                coordinates[atom_number] =\
                np.array([float(line[30:38].strip()),
                    float(line[38:46].strip()),
                    float(line[46:54].strip())])
                #print(coordinates) #bug check = ok
                #print(coordinates[atom_number])#bug check = ok
               
                '''
                res_type[atom_number] = line[17:20].strip()
                atom_type[atom_number] = atom_type
                chains[atom_number] = line[21:22].strip()
                res_number[atom_number] = int(line[22:26].strip())
                occupency[atom_number] = float(line[54:60].strip())
                temperature_factors[atom_number] = float(line[60:66].strip())
                element_symbol[atom_number] = line[76:78].strip()
                atom_charge[atom_number] = line[78:80].strip()
                '''
    print(coordinates) # bug check = youpi !    
    return(coordinates)
    

if __name__ == "__main__":                  
    parser = argparse.ArgumentParser("Manage pdp files and isolate alpha carbons")
    parser.add_argument("pdb_file", help = "file .pdb ", type=str)
    args = parser.parse_args()
    if args.pdb_file[-4:] != '.pdb' :         
        print("Please select file with pdb extension")
    else: 
        pdb_file = pdb_parser(args.pdb_file)
        #print(pdb_file) bug check