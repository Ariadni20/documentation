'''Constants'''

#AUTHOR: Ariadni Boziki

#import scipy as sc
import math
from scipy import constants

'''The constants are from scipy.constants'''

eV=constants.value('electron volt-joule relationship')
c=constants.value('speed of light in vacuum')
Ang=1.0e-10

cart_pol_factor = 0.062415091 # C/m2 to e/A2
dipole_factor = (eV/(1./(10*c)))/Ang  #(e -> D/Ang)
ir_factor = 1

bohr2m = 0.529177249e-10
angstrom2m = 1e-10
hartree2joule = 4.35974434e-18
eV2joule = 1.60218e-19
hartree2joule = 4.35974e-18
speed_of_light = 299792458
avogadro = 6.0221413e+23

# For FHI-aims
#vib_constant = math.sqrt((avogadro*eV2joule*1000)/(angstrom2m*angstrom2m))/(2*math.pi*speed_of_light*100)
# For DFTB+
vib_constant = math.sqrt((avogadro*hartree2joule*1000)/(bohr2m*bohr2m))/(2*math.pi*speed_of_light*100)

# Dictionary containing the atomic masses in amu.

masses = {
    "H": 1.00794,
    "N": 14.0067,
    "C": 12.0107,
    "O": 15.9994,
    "S": 32.065,
    "Cl": 35.453
}


# Dictionary containing all atomic Hubbard derivates (atomic units). Retrieved by https://dftb.org/parameters/download/3ob/3ob-3-1-cc.

hubbard_derivates = {
    "Br": -0.0573,
    "C": -0.1492,
    "Ca": -0.0340,
    "Cl": -0.0697,
    "F": -0.1623,
    "H": -0.1857,
    "I": -0.0433,
    "K": -0.0339,
    "Mg": -0.02,
    "N": -0.1535,
    "Na": -0.0454,
    "O": -0.1575,
    "P": -0.14,
    "S": -0.11,
    "Zn": -0.03
}

# Dictionary containing maximum angular momenta.

max_angular_mom = {
    'Br': "d",
    'C': "p",
    'Ca': "p",
    'Cl': "d",
    'F': "p",
    'H': "s",
    'I': "d",
    'K': "p",
    'Mg': "p",
    'N': "p",
    'Na': "p",
    'O': "p",
    'P': "d",
    'S': "d",
    'Zn': "d"
}
