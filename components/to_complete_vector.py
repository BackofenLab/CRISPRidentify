import numpy as np


def get_full_vector(list_vectors):
    """
    8: (2, 4, 5, 6, 7, 8, 9, 11),
    9: (1, 2, 4, 5, 7, 8, 9, 10, 12),
    10: (0, 2, 3, 4, 5, 6, 7, 10, 11, 12),
    """
    fist_model_vector = list_vectors[0][0]
    second_model_vector = list_vectors[1][0]
    third_model_vector = list_vectors[2][0]

    f_0 = third_model_vector[0]
    f_1 = second_model_vector[0]
    f_2 = fist_model_vector[0]
    f_3 = third_model_vector[2]
    f_4 = third_model_vector[3]
    f_5 = third_model_vector[4]
    f_6 = third_model_vector[5]
    f_7 = third_model_vector[6]
    f_8 = fist_model_vector[5]
    f_9 = second_model_vector[6]
    f_10 = second_model_vector[7]
    f_11 = fist_model_vector[7]
    f_12 = second_model_vector[8]

    return np.asarray([f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11, f_12]).reshape(1, -1)