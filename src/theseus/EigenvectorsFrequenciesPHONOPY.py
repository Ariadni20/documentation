'''Eigenvectors and eigenvalues from PHONOPY code'''

#AUTHOR: Ariadni Boziki

import yaml
import math
import numpy as np
import scipy.linalg as spl
from theseus import Constants as c
from theseus import InputsPreparation as inputs


def read_eig_vec_phonopy(code, dynamical_matrix):

    """It returns the eigenvalues, eigenvectors and in turn
    frequencies."""

    hessian = dynamical_matrix
    eigvals, eigvecs = np.linalg.eigh(hessian)

#   Different method for the diagonalization. The results of both methods 
#   are almost identical
#   eigvals_test, eigvecs_test = spl.eigh(hessian_mass_weighted)

    frequencies = np.sqrt(np.abs(eigvals)) * np.sign(eigvals)

#   Frequencies in 1/cm
    frequencies_in_1_over_cm = np.sqrt(np.abs(eigvals)) * np.sign(eigvals) * c.vib_constant

    if (code == 'aims'):
        conversion_factor_to_THz = 15.633302
        conversion_factor_to_cm_minus_1 = 15.633302*33.356
        frequencies_THz = frequencies * conversion_factor_to_THz
        frequencies_in_cm_minus_1 = frequencies * conversion_factor_to_cm_minus_1

    if (code == 'dftb+'):
        conversion_factor_to_THz = 154.10794
        conversion_factor_to_cm_minus_1 = 154.10794*33.356
        frequencies_THz = frequencies * conversion_factor_to_THz
        frequencies_in_cm_minus_1 = frequencies * conversion_factor_to_cm_minus_1


    return eigvecs, eigvals, frequencies_in_cm_minus_1
