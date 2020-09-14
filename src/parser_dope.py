
#!/usr/bin/env python3

import argparse
import numpy as np

def dope_parser(dope_file) : 

    '''
    This function manage dope file isolating alpha carbon - alpha carbon
    It retrieves : 
    1 - energy value for each residue - residue

    Args : 
        dope_file : Statistical potential file `dope.par` given by Jean-Christophe Gelly
        extension .txt

    Returns : 
        Dictionnary: {residue 1 - residue 2 : array[energy_value from 1 to end]}
    '''

    with open(dope_file, "r") as fillin:
        lines = fillin.readlines()
        energy_value = {}
        for line in lines : 
            word = line.split(" ")
            #print(dope_list) # check = OK
            res_first = word[0].strip()
            atom_first = word[1].strip()
            res_sec = word[2].strip()
            atom_sec = word[3].strip()
            #print(atom_first) # check = OK
            if atom_first == "CA" and atom_sec == "CA":
                #print("we are good, condition OK") # check = OK
                res_res = str(res_first + "-" + res_sec)
                #print(res_res) #check = OK
                energy_value[res_res] =\
                np.array([float(word[4].strip()),
                        float(word[5].strip()),
                        float(word[6].strip()),
                        float(word[7].strip()),
                        float(word[8].strip()),
                        float(word[9].strip()),
                        float(word[10].strip()), 
                        float(word[11].strip()),
                        float(word[12].strip()), 
                        float(word[13].strip()),
                        float(word[14].strip()),
                        float(word[15].strip()),
                        float(word[16].strip()),
                        float(word[17].strip()),
                        float(word[18].strip()),
                        float(word[19].strip()),
                        float(word[20].strip()), 
                        float(word[21].strip()),
                        float(word[22].strip()), 
                        float(word[23].strip()),
                        float(word[24].strip()),
                        float(word[25].strip()),
                        float(word[26].strip()),
                        float(word[27].strip()),
                        float(word[28].strip()),
                        float(word[29].strip()),
                        float(word[30].strip()), 
                        float(word[31].strip()),
                        float(word[32].strip()), 
                        float(word[33].strip())])
        print(energy_value)
        return(energy_value)


if __name__ == "__main__":                  
    parser = argparse.ArgumentParser("Manage the dope.par.txt file : isolate CA-CA and energy values")
    parser.add_argument("dope_file", help = "dope.par.txt ", type=str)
    args = parser.parse_args()
    if args.dope_file[:] != 'dope.par.txt' :         
        print("Please select only the file 'dope.par.txt' provided by Jean-Christophe Gelly")
    else: 
        dope_file = dope_parser(args.dope_file)


