'''It converts the coordinates from cartesian to fractional.'''

#AUTHOR: Ariadni Boziki

import numpy as np
from theseus import InputsPreparation as inputs
import os
from theseus import FiniteDisplacements_phonopy as disps


def vector_inside_the_box(super_invmat, m, super_mat):

    '''It moves the atoms into the simulation box.'''


    real_vectors = np.dot(super_invmat,m)
    vectors_in_the_box = (real_vectors + 1e-5) % 1.0 - 1e-5
    dot_vectors_in_the_box = np.dot(super_mat,vectors_in_the_box)

    return dot_vectors_in_the_box


def cartesian_to_fractional_DFTBplus(supercell=False):

    '''It converts cartesian coordinates to fractional - DFTB+.'''

    if supercell:
        geo_input = 'geo.genS'
    else:
        geo_input = 'geo.gen'

    geometry_processor = inputs.GeometryProcessor(geo_input, 'dftb+')

    lattice = geometry_processor.read_lattice()
    positions = geometry_processor.read_coordinates()
    numbers, types = geometry_processor.read_the_atom_type()
    no_atoms = geometry_processor.number_of_atoms() 

    geometry_input = open(geo_input, 'r')
    atom_types_dftb = geometry_input.readlines()[1]
    atom_types_dftb = atom_types_dftb.split()

    inputs.copy_files_in_the_same_directory('geo.gen', 'geo_backup.gen')

    axx = lattice[0][0]
    axy = lattice[0][1]
    axz = lattice[0][2]
    bxx = lattice[1][0]
    bxy = lattice[1][1]
    bxz = lattice[1][2]
    cxx = lattice[2][0]
    cxy = lattice[2][1]
    cxz = lattice[2][2]

    lattice = np.array(lattice)
    lattice_inv_mat = np.linalg.inv(lattice)
    super_mat = lattice.transpose()
    super_invmat = np.linalg.inv(super_mat)

    path = os.getcwd()
    output = os.path.join(path, "geo.gen.F")
    with open(output, 'w') as ff:
        ff.writelines('%s F\n' %(no_atoms))
        for i in atom_types_dftb:
            ff.writelines('%s ' '' %(i))
        ff.writelines('\n')
        for i, j, k in zip(positions, numbers, range(no_atoms)):

            kk = k + 1

            ii = np.array(i)
            trans_conv = vector_inside_the_box(super_invmat, ii, super_mat)
            trans_conv = np.dot(super_invmat,trans_conv)

            ff.writelines(' ' ' %s ' ' %s ' ' %s ' ' %s ' ' %s \n' %(kk, j, trans_conv[0], trans_conv[1], trans_conv[2]))
        ff.writelines(' ' ' 0.0 ' ' 0.0 ' ' 0.0\n')
        ff.writelines(' ' ' %s ' ' %s ' ' %s\n' %(axx, axy, axz))
        ff.writelines(' ' ' %s ' ' %s ' ' %s\n' %(bxx, bxy, bxz))
        ff.writelines(' ' ' %s ' ' %s ' ' %s\n' %(cxx, cxy, cxz))

    disps.rename('geo.gen.F', 'geo.gen')


def cartesian_to_fractional_FHIaims(supercell=False):


    '''It converts cartesian coordinates to fractional - FHI-aims.'''

    if supercell:
        geo_input = 'geometry.in.supercell'
    else:
        geo_input = 'geometry.in'

    geometry_processor = inputs.GeometryProcessor(geo_input, 'aims')

    lattice = geometry_processor.read_lattice()
    positions = geometry_processor.read_coordinates()
    numbers, types = geometry_processor.read_the_atom_type()
    no_atoms = geometry_processor.number_of_atoms() 

    inputs.copy_files_in_the_same_directory('geometry.in', 'geometry_backup.in')

    axx = lattice[0][0]
    axy = lattice[0][1]
    axz = lattice[0][2]
    bxx = lattice[1][0]
    bxy = lattice[1][1]
    bxz = lattice[1][2]
    cxx = lattice[2][0]
    cxy = lattice[2][1]
    cxz = lattice[2][2]

    lattice = np.array(lattice)
    lattice_inv_mat = np.linalg.inv(lattice)
    super_mat = lattice.transpose()
    super_invmat = np.linalg.inv(super_mat)

    path = os.getcwd()
    output = os.path.join(path, "geometry.in.F")
    with open(output, 'w') as ff:
        ff.writelines('\n')
        ff.writelines(' lattice_vector   %s  %s  %s\n' %(axx, axy, axz))
        ff.writelines(' lattice_vector   %s  %s  %s\n' %(bxx, bxy, bxz))
        ff.writelines(' lattice_vector   %s  %s  %s\n' %(cxx, cxy, cxz))
        ff.writelines('\n')
        for i, j in zip(positions, types):

            ii = np.array(i)
            trans_conv = vector_inside_the_box(super_invmat, ii, super_mat)
            trans_conv = np.dot(super_invmat,trans_conv)

            ff.writelines(' atom_frac  %s  %s  %s  %s\n' %(trans_conv[0], trans_conv[1], trans_conv[2], j))

    disps.rename('geometry.in.F', 'geometry.in')
