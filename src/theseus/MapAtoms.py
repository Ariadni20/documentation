'''It takes care of equivalent atoms. It maps the polarizability and 
cartesian polarization of irreducible atoms back to reducible atoms.''' 

#AUTHOR: Ariadni Boziki

import spglib
import numpy as np
from theseus import InputsPreparation as inputs
from theseus import Coordinates as coords


def equivalent_symmetry_atoms(code):

    """It returns a dictionary with the equivalent atoms.
    Spglib is used for this reason. This is equivalent to what
    PHONOPY does."""

    equiv_atoms = dict()	

    if code == 'aims':
        geometry_processor = inputs.GeometryProcessor('geometry.in', 'aims')
    if code == 'dftb+':
        geometry_processor = inputs.GeometryProcessor('geo.gen', 'dftb+')
        
    lattice = geometry_processor.read_lattice()
    positions = geometry_processor.read_coordinates()
    numbers, types = geometry_processor.read_the_atom_type()

    cell = (lattice, positions, numbers)

    dataset = spglib.get_symmetry_dataset(cell)

    equivalent_atoms = dataset['equivalent_atoms']

    index = 0 
    for item in equivalent_atoms:
        equiv_atoms[index] = item
        index += 1	

    return equivalent_atoms


def map_axis_coord(coord):

	"""It returns the coordinate, (x, y, z)
	for mapping the equivalent coordinates.
	"""

	axis_coord = 0	

	if (coord=='x'):
		axis_coord=1
	if (coord=='y'):
		axis_coord=2
	if (coord=='z'):
		axis_coord=3

	return axis_coord


def merge_arrays(pol, cartesian_pol, element_axis_coord):

	"""It adds in the polarizablty and cartesian polarization
	arrays two more columns. One with the number of the atom
	and a second one with the coordinate."""

	element_axis_coord_float = element_axis_coord.astype(np.float)

	pol_extended = np.concatenate((pol, element_axis_coord_float), axis=1)
	cartesian_pol_extended = np.concatenate((cartesian_pol, element_axis_coord_float), axis=1)

	return pol_extended, cartesian_pol_extended


def map_properties_of_equivalent_atoms(pol_extended, cartesian_pol_extended, code):

    """It returns polarizability and cartesian polarization arrays
    after having mapped the polarizabilty and cartesian polarization 
    of reducible atoms with irreducible atoms (equivalent atoms).
    It also considers the case where only one coordinate of an atom
    has been displaced."""

    pol = np.empty([0,8])
    cartesian_pol = np.empty([0,5])

    equivalent_atoms = equivalent_symmetry_atoms(code)

    for atom, equiv_atom in enumerate(equivalent_atoms):
        for values in pol_extended:
            if (values[6] == equiv_atom):
                if (values[7] == 1.0):
                    pol = np.append(pol,[values], axis = 0)
        for values in pol_extended:
            if (values[6] == equiv_atom):
                if (values[7] == 2.0):
                    pol = np.append(pol,[values], axis = 0)
        for values in pol_extended:
            if (values[6] == equiv_atom):
                if (values[7] == 3.0):
                    pol = np.append(pol,[values], axis = 0)


        for values in cartesian_pol_extended:
            if (values[3] == equiv_atom):
                if (values[4] == 1.0):
                    cartesian_pol = np.append(cartesian_pol,[values], axis = 0)
        for values in cartesian_pol_extended:
            if (values[3] == equiv_atom):
                if (values[4] == 2.0):
                    cartesian_pol = np.append(cartesian_pol,[values], axis = 0)
        for values in cartesian_pol_extended:
            if (values[3] == equiv_atom):
                if (values[4] == 3.0):
                    cartesian_pol = np.append(cartesian_pol,[values], axis = 0)


    length_pol = len(pol)

    if code == 'aims':
        geometry_processor = inputs.GeometryProcessor('geometry.in', 'aims')
    if code == 'dftb+':
        geometry_processor = inputs.GeometryProcessor('geo.gen', 'dftb+')

    number_of_atoms = inputs.number_of_atoms()

    if number_of_atoms != length_pol:
        for ii in range(number_of_atoms):
            count = (pol[:,6] == ii).sum()
            if count != 0:
                if(count % 3 != 0):
                    check = inputs.check_similarities_of_coord(geometry_processor.read_coordinates(),ii)
                    if check:
                        if count == 1:
                            kk = 0
                            for jj in pol:
                                if jj[6] == ii:
                                    line = kk	
                                kk += 1
                            pol = np.insert(pol, line, pol[line,:], 0)
                            pol[line+1, 7] = 2	
                            pol = np.insert(pol, line+1, pol[line+1,:], 0)
                            pol[line+2, 7] = 3	

        for ii in range(number_of_atoms):
            count = (cartesian_pol[:,3] == ii).sum()
            if count != 0:
                if(count % 3 != 0):
                    check = inputs.check_similarities_of_coord(geometry_processor.read_coordinates(),ii)
                    if check:
                        if count == 1:
                            kk = 0
                            for jj in cartesian_pol:
                                if jj[3] == ii:
                                    line = kk	
                                kk += 1
                            cartesian_pol = np.insert(cartesian_pol, line, cartesian_pol[line,:], 0)
                            cartesian_pol[line+1, 4] = 2	
                            cartesian_pol = np.insert(cartesian_pol, line+1, cartesian_pol[line+1,:], 0)
                            cartesian_pol[line+2, 4] = 3	


        return pol, cartesian_pol


def remove_columns_from_array(pol_extended, cartesian_pol_extended):

	"""It returns the polarizability and cartesian polarization arrays
	after having removed the last two columns that correspond to the
	number and coordinate (axis) of the atom.
	"""

	pol = np.delete(pol_extended,np.s_[6:8],axis=1)
	cartesian_pol = np.delete(cartesian_pol_extended,np.s_[3:5],axis=1)


	return pol, cartesian_pol


def get_international_space_group_number(code):


    """It returns the international space group number"""

   
    if (code == 'aims'):
        geometry_processor = inputs.GeometryProcessor('geometry.in', 'aims')
        with open('geometry.in','r') as file:
            file_content = file.read()
            if 'atom_frac' in file_content:
                pass
            else:
                coords.cartesian_to_fractional_FHIaims(supercell=False)
                
    if (code == 'dftb+'):
        geometry_processor = inputs.GeometryProcessor('geo.gen', 'dftb+')
        with open('geo.gen','r') as file:
            file_content = file.read()
            if 'F' in file_content:
                pass
            else:
                coords.cartesian_to_fractional_DFTBplus(supercell=False)
            
    lattice = geometry_processor.read_lattice()
    positions = geometry_processor.read_coordinates()
    numbers, types = geometry_processor.read_the_atom_type()

    cell = (lattice, positions, numbers)

    dataset = spglib.get_symmetry_dataset(cell)

    space_group_number = dataset['number']

    return space_group_number
