# DOPE energy values spread from 0.25 to 15.25A by 0.5 intervals (30 intervals)



import numpy as np
distance = np.array([0.45, 2.89, 5.9, 0.56, 4.10])
print(distance[1])

dope_dict = \
{('ALA', 'ALA'): [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, -1.64, 1.02, -0.0, -1.23,\
                  -0.62, -0.63, 0.33, 0.58, 0.48, 0.1, -0.33, -0.07, -0.3, -0.25, -0.16, 0.02,\
                  0.03, -0.08, 0.01, -0.02, -0.08, -0.12, -0.02]}

interval_index = round(int(distance[1] * 2 - 0.5))
print(dope_dict['ALA', 'ALA'][interval_index])
energy[i, j] = dope_dict[('residu_fixe_i', 'les autres_j')][interval_index]



