'''IR and Raman intensities. For periodic systems.'''

#AUTHOR: Ariadni Boziki

import numpy as np
from theseus import Constants as c
from theseus import InputsPreparation as inputs

def IR_intensity(volume, cartesian_pol, eig_vec):

	"""It returns the Infrared Intensity.""" 

	grad_cart_pol = c.cart_pol_factor * c.dipole_factor * volume * cartesian_pol # D/Ang
	IR_intensity = np.sum(np.dot(np.transpose(grad_cart_pol), eig_vec)**2, axis=0) * c.ir_factor # D^2/A^2*amu

#	print(f'IR_intensity: {IR_intensity}')

	return IR_intensity


def Raman_activity(pol, eig_vec, code):

    """It returns the Raman Activity."""

    alphas = np.dot(np.transpose(pol), eig_vec)

    if (code == 'aims'):
        geometry_processor = inputs.GeometryProcessor('geometry_backup.in', 'aims') 
    if (code == 'dftb+'):
        geometry_processor = inputs.GeometryProcessor('geo_backup.gen', 'dftb+') 
    
    natoms = geometry_processor.number_of_atoms()
   
    xx = np.zeros(natoms*3)
    yy = np.zeros(natoms*3)
    zz = np.zeros(natoms*3)
    xy = np.zeros(natoms*3)
    xz = np.zeros(natoms*3)
    yz = np.zeros(natoms*3)

    for step in range(natoms*3):	
        xx[step] = alphas[0,step]
        yy[step] = alphas[1,step]
        zz[step] = alphas[2,step]
        xy[step] = alphas[3,step]
        xz[step] = alphas[4,step]
        yz[step] = alphas[5,step]

    alpha = (xx + yy + zz) * (1./3)
    beta = (xx - yy)**2 + (xx - zz)**2 + (yy - zz)**2 + 6*(xy**2 + xz**2 + yz**2)

    raman_activity = 45 * (alpha**2) + (7./2) * beta
    raman_activity = raman_activity * 0.02195865620442408 #bohr^6/ang^2 to ang^4

#    print(f'Raman activity: {raman_activity}')

    return raman_activity
