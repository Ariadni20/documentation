'''Preparation and modification of input files'''

#AUTHOR: Ariadni Boziki

import os
import shutil
import subprocess
import glob

class GeometryProcessor():

    """
    Class for exctracting information out of a geometry input file.

    Attributes
    ----------

    geom_input: str
        Namefile of the geometry input file.
    code: str
        Code
    """

    def __init__(self, geom_input: str, code: str):

        self.geom_input = geom_input
        self.code = code


    def number_of_atoms(self):

        """
        It returns the number of atoms in the structure.

        index: int
            Nunmber of atoms in a given geometry input file.
        """

        index = 0
        with open(self.geom_input,'r') as geometry_input:
            for lines in geometry_input:
                if self.code == 'aims' and any(keyword in lines for keyword in ('atom', 'atom_frac')):
                    index += 1
                elif self.code == 'dftb+':
                    no_of_columns = len(lines.split())
                    if no_of_columns == 5:
                        index += 1
        return index


    def read_lattice(self):

        """
        It returns the lattice cell (9 elements) from the geometry input file.

        cell: list
            List containing the lattice cell elements as floats.
        """
	
        cell = []

        with open(self.geom_input, 'r') as geometry_input:

            if self.code == 'aims':
                for lines in geometry_input:
                    if 'lattice_vector' in lines:
                        cell_tmp = lines.split()[1:]
                        cell_tmp_float = [float(s.replace(',','')) for s in cell_tmp]
                        cell.append(cell_tmp_float)

            elif self.code == 'dftb+':
                for line in (geometry_input.readlines() [-3:]):
                    cell_tmp = line.split()
                    cell_tmp_float = [float(s.replace(',','')) for s in cell_tmp]
                    cell.append(cell_tmp_float)

        return cell


    def read_coordinates(self):
    
        """
        It returns the coordinates from the geometry input file.
        
        coord : list
            List containing coordinate tuples as floats.
        """

        coord = []	

        with open(self.geom_input, 'r') as geometry_input:

            if self.code == 'aims':
                for lines in geometry_input:
                    if any(keyword in lines for keyword in ('atom', 'atom_frac')):
                        coord_tmp = lines.split()[1:4]
                        coord_tmp_float = [float(s.replace(',','')) for s in coord_tmp]
                        coord.append(coord_tmp_float)

            elif self.code == 'dftb+':
                for lines in geometry_input:
                    no_of_columns = len(lines.split())
                    if no_of_columns == 5:
                        coord_tmp = lines.split()[2:5]
                        coord_tmp_float = [float(s.replace(',','')) for s in coord_tmp]
                        coord.append(coord_tmp_float)

        return coord


    def read_the_atom_type(self):

        """
        It returns the atom type of the atoms of the structure from the geometry input file.

        atom_type_no: tuple
            Tuple containing lists of atom type numbers.
        atom_type: tuple
            Tuple containing lists of atom types.
        """

        atom_type = []
        atom_type_no = []
        atom_types_tmp_dftb = []

        with open(self.geom_input,'r') as geometry_input:

            contents = geometry_input.readlines()

            if self.code == 'dftb+':
                atom_types_tmp_dftb = contents[1].split()
                
            for lines in contents:
                
                if self.code == 'aims' and any(keyword in lines for keyword in ('atom', 'atom_frac')):
                    atom_type_tmp = lines.split()[-1]
                    atom_type.append(atom_type_tmp)
                    
                    if len(atom_type) == 1:
                        atom_type_no_tmp = 1
                    elif atom_type[-1] == atom_type[-2]:
                        atom_type_no_tmp = atom_type_no[-1]
                    else:
                        atom_type_no_tmp = atom_type_no[-1] + 1
                    
                    atom_type_no.append(atom_type_no_tmp)
                
                elif self.code == 'dftb+':
                    no_of_columns = len(lines.split())
                    if no_of_columns == 5:
                        atom_type_no_tmp = lines.split()[1]
                        atom_type_no.append(atom_type_no_tmp)
            if (self.code == 'dftb+'):
                for ii in range(len(atom_types_tmp_dftb)):
                    kk = ii+1
                    for jj in atom_type_no:
                        ll = int(jj)
                        if kk is ll:
                            atom_type.append(atom_types_tmp_dftb[ii])

        return atom_type_no, atom_type



def copy_files_in_dir(pattern, filename):
	
    """Copy input files in the generated directories."""

    contents = os.listdir(os.getcwd())
    for item in contents:
        if os.path.isdir(item):
            shutil.copy(filename, item)	


def copy_files_in_the_same_directory(filename, new_name):

    '''It copies a file in the same directory.'''

    shutil.copy(filename, new_name)


def check_similarities_of_coord(coord,no_of_line):

	"""It returns True if the coordinates are the same."""

	line = coord[no_of_line]
	if (line[0] == line[1]) and ((line[0] == line[2])) and (line[1] == line[2]):
		return True
	else:
		return False


def _check_similarities_of_coord(coord, no_of_line):

    """
    It checks if the coordinates are the same.

    Returns:
    book: True if the coordinates are the same, False otherwise.
    """

    line = coord[no_of_line]
    return all(x == line[0] for x in line[1:]) 
