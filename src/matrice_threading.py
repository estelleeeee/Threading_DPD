import math
import numpy as np 

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
