
#!/usr/bin/env python3

import argparse

def dope_parser(dope_file) : 

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
        for line in fillin : 
            words_list = line.strip().split(" ")
            res_first = words_list[0].strip()
            atom_first = words_list[1].strip()
            res_sec = words_list[2].strip()
            atom_sec = words_list[3].strip()
            if atom_first == "CA" and atom_sec == "CA":
                res_res = (res_first, res_sec)
                bin_list = [float(i) for i in words_list[4:]]
                energy_value_dict[res_res] = bin_list
        #print(energy_value_dict[("ALA","ALA")])
        return(energy_value_dict)


if __name__ == "__main__":                  
    parser = argparse.ArgumentParser("Manage the dope.par.txt file to isolate CA-CA pair energies")
    parser.add_argument("-d", "--dope-file", help = "dope.par.txt ", type=str)
    args = parser.parse_args()
    DOPE_FILE_NAME = args.pdb_file
    if DOPE_FILE_NAME[:] != 'dope.par.txt' :         
        print("Please select only the file 'dope.par.txt' provided by Jean-Christophe Gelly")
    else: 
        dope_energy_value = dope_parser(DOPE_FILE_NAME)


