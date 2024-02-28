'''Two-point central difference to calculate the polarizability
and cartesian polarization. For periodic systems.
'''

#AUTHOR: Ariadni Boziki

import os
import re
import itertools as it
import numpy as np
from theseus import MapAtoms as mapatoms


def subgroups_finite_disp_drct():

	'''It returns a list whose rows consist of the names of two geometry files.
	In these two files the same coordinate has been displaced.
	'''
	
	drct = []
	contents = filter(os.path.isdir, os.listdir(os.getcwd()))
	for i in contents:
		if 'Coord' in i:
			drct.append(i)
			
	groupsdir = []
	
	for dirname, item in it.combinations(drct, 2):
		
		atom_coord1 = dirname.split('-')[1]
		atom_coord2 = item.split('-')[1]

		if atom_coord1 == atom_coord2:
			save = dirname, item
			groupsdir.append(save)

	return groupsdir


def split_line(line):

	'''It splits the lines and saves the values in numpy arrays.'''

	larray=np.array(line.strip().split(' '))
	values=larray[larray!='']
	return values


def find_pattern(output_file):

	'''It finds specific patterns in FHI-aims output files
	and returns the values. Used for Cartesian polarization 
	and Polarizability.'''

	myfile = open(output_file, 'r')
	for line in myfile:
		if 'Cartesian Polarization' in line:
			CartesianPolarization = np.float64(split_line(line)[-3:])
		if 'Polarizability' in line:
			Polarizability = np.float64(split_line(line)[-6:])

	return CartesianPolarization, Polarizability


def find_pattern_DFTBplus(output_file, dispersion=False):

    '''It finds specific patterns in DFTB+ output files
    and returns the values. Used for Dipole moment
    and Polarizability.'''

    path = os.getcwd()

    if dispersion:
        drct = 'polarizability'
        path_drct = os.path.join(path, drct)
        os.chdir(path_drct)
    else:
        os.chdir(path)
    with open(output_file, 'r') as f:
        lines = f.readlines()
        for index, line in enumerate(lines):
            if 'Static polarisability:' in line:
                pol1 = lines[index+1]
                pol1 = np.float64(split_line(pol1))
                pol2 = lines[index+2]
                pol2 = np.float64(split_line(pol2))
                pol3 = lines[index+3]
                pol3 = np.float64(split_line(pol3))

                xx = pol1[0]
                yy = pol2[1]
                zz = pol3[2]
                xy = pol1[1]
                xz = pol1[2]
                yz = pol2[2]

                Polarizability = np.array([xx, yy, zz, xy, xz, yz])

    ff =  open('detailed.out', 'r')
    for line in ff:
        if 'Dipole moment:' in line:
            CartesianPolarization = np.float64(split_line(line)[-4:-1])

    return CartesianPolarization, Polarizability


def UnitCellVolume(output_file):

    '''It returns the unit cell volume. Given that during the frozen 
    phonon approximation, no unit cell and geometry optimizations are applied
    the volume can be extracted by any FHI-aims output.'''

    myfile = open(output_file, 'r')
    for line in myfile:
        if '| Unit cell volume' in line:
            UnitCellVolume = np.float64(split_line(line)[-2])
        if 'Volume:' in line:
            UnitCellVolume = np.float64(split_line(line)[-2])

    return UnitCellVolume


def pol_fhiaims(sign_atom_coord, disp, code, dispersion, output_file):

    '''Two-point central difference is applied as the numerical differentiation
    method. The two structures in which the same coordinate has been displaced
    are used in each step.
    '''

    path = os.getcwd()
    c_diff_fraction = - 1. / (2. * disp)

    pol = np.empty([0,6])	
    cartesian_pol = np.empty([0,3])
    element_axis_coord = np.empty([0,2])

    for item in sign_atom_coord:
        dir1 = item[0]      # First directory
        sign1 = dir1[-1]
        no_of_element = dir1.split('-')[2]
        no_axis_coord = mapatoms.map_axis_coord(dir1.split('-')[3]) 
        path_dir1 = os.path.join(path, dir1)
        os.chdir(path_dir1)
        if (code == 'aims'):
            cart_pol1, polarizability1 = find_pattern(output_file)
        if (code == 'dftb+'):
            cart_pol1, polarizability1 = find_pattern_DFTBplus(output_file, dispersion)
        if dispersion:
            os.chdir('../../')
        else:
            os.chdir('../')

        dir2 = item[1]     # Second directory
        sign2 = dir2[-1]
        path_dir2 = os.path.join(path, dir2)
        os.chdir(path_dir2)
        if (code == 'aims'):
            cart_pol2, polarizability2 = find_pattern(output_file)
        if (code == 'dftb+'):
            cart_pol2, polarizability2 = find_pattern_DFTBplus(output_file, dispersion)
        if dispersion and (code == 'dftb+'):
            os.chdir('../../')
        else:
            os.chdir('../')

        if sign1 == '+':
            coeff1 = 1
        else:
            coeff1 = -1

        if sign2 == '+':
            coeff2 = 1
        else:
            coeff2 = -1

#       Polarizability
        pol_tmp = polarizability1 * coeff1 * c_diff_fraction + polarizability2 * coeff2 * c_diff_fraction
        pol = np.append(pol,[pol_tmp], axis = 0)
#       Cartesian Polarization
        cartesian_pol_tmp = cart_pol1 * coeff1 * c_diff_fraction + cart_pol2 * coeff2 * c_diff_fraction
        cartesian_pol = np.append(cartesian_pol,[cartesian_pol_tmp], axis = 0)

        tmp_element_axis_coord = [no_of_element, no_axis_coord]
        element_axis_coord = np.append(element_axis_coord, [tmp_element_axis_coord], axis = 0) 

    return pol, cartesian_pol, element_axis_coord
