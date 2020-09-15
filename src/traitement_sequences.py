import numpy as np 

def code_1_lettre_3_lettres(fasta_file):
    """Fonction permettant de transformer un fichier fasta avec AA en code 1 lettre en 
    un dictionnaire avec AA, code 3 lettres, sous forme d'une liste.
    
    Parameters
    ----------
    fasta_file : fichier fasta
    
    Returns
    -------
    dict[str]=[str,str,...]
    """
    DICO_CONVERSION_AA = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', \
    'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', \
    'P':'PRO', 'S':'SER','T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
    dico_AA_3lettres = {}
    with open(fasta_file, "r") as file:
        seq_AA_3lettres = []
        for line in file:
            if line[0] == ">":
                nom_seq = line[1:5] #le nom de la sequence sera notre cle -> code en dur
            else:
                for AA in line[:-1]:
                    if AA == 'X': #HETATOM (ne nous interesse pas)
                        pass
                    else:
                        seq_AA_3lettres.append(DICO_CONVERSION_AA[AA])
    dico_AA_3lettres[nom_seq] = seq_AA_3lettres
    return dico_AA_3lettres

def pdb_parser(pdb_file):
    """This function manage pdb file isolating alpha carbons only (atome type).
    It retrieves:
    1. residu number
    2. x, y, z coordinates
    
    Args
    ----
        pdb_file : file with pdb extension

    Returns
    -------
        Dictionnary: {atom_number : array[coord_x, coord_y, coord_z]}
    """

    with open(pdb_file, "r") as fillin:
        coordinates_dict = {}
        for line in fillin :
            atom_type = line[12:16].strip()
            if line[0:6].strip() == "ATOM" and (atom_type == "CA" or atom_type == "C1*"):
                res_number = int(line[22:26].strip())
                coordinates_dict[res_number] =\
                    np.array([float(line[30:38].strip()),
                              float(line[38:46].strip()),
                              float(line[46:54].strip())])
    return(coordinates_dict)
