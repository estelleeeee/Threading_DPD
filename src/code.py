import math
import numpy as np

def code_1_lettre_3_lettres(fasta_file):
    """Fonction permettant de transformer un fichier fasta avec AA en code 1 lettre en 
    un dictionnaire avec AA, code 3 lettres, sous forme d'une liste.
    
    Parameters
    ----------
    fasta_file : fichier fasta
    
    Returns
    -------
    liste de str code 3 lettres
    """
    DICO_CONVERSION_AA = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', \
    'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', \
    'P':'PRO', 'S':'SER','T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
    list_AA_3lettres = []
    with open(fasta_file, "r") as file:
        for line in file:
            if line[0] == ">":
                nom_seq = line[1:5] #le nom de la sequence sera notre cle -> code en dur
            else:
                for AA in line[:-1]:
                    if AA == 'X': #HETATOM (ne nous interesse pas)
                        pass
                    else:
                        list_AA_3lettres.append(DICO_CONVERSION_AA[AA])
    return list_AA_3lettres

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

def matrice_distance(dico_AA):
    """Creer une matrice de distance entre tous les AA du dictionnaire

    Parameters
    ----------
    dico_AA : dict
        contient la position et les coordonnees d'un AA

    Returns
    -------
    np.array(matrice)
        matrice contenant les distances entre AA
    """
    matrice = []
    for num1 in dico_AA:
        dist_ligne = []
        for num2 in dico_AA:
            dist_ligne.append(distance_array(dico_AA[num1],dico_AA[num2]))
        matrice.append(dist_ligne)
    return np.array(matrice)

def distance_array(array1,array2):
    """Calcul la distance entre deux arrays : [x, y, z]
    
    Parameters
    ----------
    array1 : numpy.ndarray
        coordonnees de la premiere position
    array2 : numpy.ndarray
        coordonnees de la deuxieme position
    
    Returns
    -------
    np.array(dist)
        distance en array entre les 2 arrays
    """
    x = array1[0] - array2[0]
    y = array1[1] - array2[1]
    z = array1[2] - array2[2]
    dist = math.sqrt(x**2 + y**2 + z**2)
    return np.array(dist)

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

def dope_score(res1, res2, dist_res1_res2):
    dope_dict = dope_parser('dope.par.txt')
    interval_index = round(int(dist_res1_res2 * 2 - 0.5))
    return dope_dict[(res1, res2)][interval_index]

def mat_dope_res_pos(dist_mat, query, res_i, pos):
    """Realise la matrice de bas niveau pour une position et un AA choisis

    Parameters
    ----------
    mat_dist: numpy.ndarray
        matrice de distance entre les coordonnees du modele
    list_seq_cible: dict
        dictionnaire contenant la sequence cible
    AA_choisi: str
        l'AA que nous allons fixer
    position: int
        la position ou nous allons fixer l'AA
    
    Returns
    -------
    
    dist_pos = dist_mat[pos-1] #distance entre AA et les positions
    dope_mat = np.empty((len(query)+1,len(dist_pos)+1))
    dope_mat.fill(np.nan)
    dope_mat[0,:] = 0
    dope_mat[:,0] = 0
    for num_AA in range(len(query)):
        if num_AA == pos: #on ignore l'AA si c'est celui fixe
            pass
        else:
            for num_dist in range(len(dist_pos)):
                if dist_pos[num_dist] == 0: #position de fixation de l'AA
                    pass 
                else:
                    # print('numA',num_AA, 'distance AA',distAA_AA)
                    dope_mat[num_AA+1][num_dist+1] = dope_score(res_i,query[num_AA],dist_pos[num_dist])
    return dope_mat
    """
    dist_pos = dist_mat[pos-1] #distance entre AA et les positions
    dope_mat = np.empty((len(query)+1,len(dist_pos)+1))
    dope_mat.fill(np.nan)
    dope_mat[0,:] = 0
    dope_mat[:,0] = 0
    #gap = 0
    for row_prev in range(1,pos):
        for col_prev in range(1,pos):
            dope_mat[row_prev][col_prev] = min(\
            dope_mat[row_prev-1][col_prev], \
            dope_mat[row_prev][col_prev-1], \
            dope_mat[row_prev-1][col_prev-1]+dope_score(res_i, query[row_prev-1], dist_pos[col_prev-1]))
    dope_mat[pos][pos] = dope_mat[pos-1][pos-1] #AA fixe, score de match = 0
    for row_after in range(pos+1, len(query)+1):
        for col_after in range(pos+1, len(dist_pos)+1):
            if math.isnan(dope_mat[row_after-1][col_after]) and math.isnan(dope_mat[row_after][col_after-1]):
                dope_mat[row_after][col_after] = dope_mat[row_after-1][col_after-1] + \
                                                dope_score(res_i, query[row_after], dist_pos[col_after])
            elif math.isnan(dope_mat[row_after-1][col_after-1]) and math.isnan(dope_mat[row_after-1][col_after]):
                dope_mat[row_after][col_after] = dope_mat[row_after][col_after-1]
            elif math.isnan(dope_mat[row_after-1][col_after-1]) and math.isnan(dope_mat[row_after][col_after-1]):
                dope_mat[row_after][col_after] = dope_mat[row_after-1][col_after]
            else: 
                dope_mat[row_after][col_after] = min(\
            dope_mat[row_after-1][col_after], \
            dope_mat[row_after][col_after-1], \
            dope_mat[row_after-1][col_after-1]+dope_score(res_i, query[row_after-1], dist_pos[col_after-1]))
    print("For the fixed residue is {:s} in position {:d}, the optimized score is {:.2f}".format(res_i, pos, dope_mat[-1][-1]))
    return dope_mat

    


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
    
