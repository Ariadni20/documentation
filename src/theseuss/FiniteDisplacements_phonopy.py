"""Creation of the subdirectories of the finite displaced geometries."""

#AUTHOR: Ariadni Boziki

import os
import shutil
import glob
from pathlib import Path
import re
import difflib
import itertools as it
import math


def even_odd_numbers(num):

	"""It returns the sign of the displacement."""

	if (num % 2) == 0:
		sign = '-'
	else:
		sign = '+'
	return sign


def diff_lines_threshold(file1, file2, code):

    """It compares two lines with three records each. In case   
    the difference between two records is 0 and the difference
    between the third record is either 0.02 for aims and 0.01 
    for dftb+ respectively, the function returns the axis of the
    displeced coordinate, (x or y or z).
    """

    axis = ''


    if (code == 'aims'):
        line1 = file1.split()[1:4]
    if (code == 'dftb+'):
        line1 = file1.split()[2:5]
    list1_tmp = [float(i) for i in line1]
    list1 = [float("{0:0.3f}".format(i)) for i in list1_tmp]
    if (code == 'aims'):
        line2 = file2.split()[1:4]
    if (code == 'dftb+'):
        line2 = file2.split()[2:5]
    list2_tmp = [float(i) for i in line2]
    list2 = [float("{0:0.3f}".format(i)) for i in list2_tmp]

    diff_1_tmp = list1[0] - list2[0]
    diff_2_tmp = list1[1] - list2[1]
    diff_3_tmp = list1[2] - list2[2]

    diff_1 = round(diff_1_tmp, 2)
    diff_2 = round(diff_2_tmp, 2)
    diff_3 = round(diff_3_tmp, 2)

    diff_1 = abs(diff_1)
    diff_2 = abs(diff_2)
    diff_3 = abs(diff_3)

    if (code == 'aims'):
        value_disp = 0.02
    if (code == 'dftb+'):
        value_disp = 0.01

    if (diff_1 == value_disp):
        axis = 'x'

    if (diff_2 == value_disp):
        axis = 'y'

    if (diff_3 == value_disp):
        axis = 'z'

    return axis


def diff_lines(file1, file2, code):

    """It compares two files line by line and returns the line
    as a number if the two lines are different. 
    It also returns the axis of the displacement.
    This information is used later on for mapping the equivalent atoms.
    """

    axis = ''
    no_of_line = 0
    final_line_1 = ''
    final_line_2 = ''
    final_no_of_line = 0

    for line_file1, line_file2 in zip(file1, file2):
        if line_file1 != line_file2:
            final_line_1 = line_file1
            final_line_2 = line_file2
            axis = diff_lines_threshold(final_line_1, final_line_2, code)
            final_no_of_line = no_of_line
        no_of_line += 1
    return  final_no_of_line, axis


def displaced_element(no_of_line, code):

    """In FHIaims input file - geometry.in, the first lines
    correspond to the cell. The atoms - coordinates start appearing
    from the 6th line.
    In DFTB+ input file - geo.gen, the first lines correspond to a comment
    and number of atoms. The atoms - coordinates start appearing
    from the 3th line.
    """

    if (code == 'aims'):
        number = 5
    if (code == 'dftb+'):
        number = 2
        
    no_of_element = no_of_line - number

    return no_of_element 


def displacements(input_file_name, code):

    """It returns a dictionary with the geometry files with the same 
    displaced coordinate. It also returns the number and coordinate of atom
    that has been displaced.
    """
	
    case_dict = dict()
    index = 0
    contents = os.listdir(os.getcwd())
    contents_list = []
    for i in contents:
        if input_file_name in i:
            contents_list.append(i)
    contents_list = sorted(contents_list)
    chunk = [contents_list[x:x+2] for x in range(0, len(contents_list), 2)]

    for item in chunk:
        file1 = open(item[0], 'r')
        file2 = open(item[1], 'r')
        no_of_line, axis = diff_lines(file1, file2, code)
        if axis == '':
            continue
        no_of_element = displaced_element(no_of_line, code)
        no1 = re.sub('[^0-9]', '', item[0])
        no2 = re.sub('[^0-9]', '', item[1])
        no1 = int(no1)
        no2 = int(no2)
        case_dict[index] = no1, no2, no_of_element, axis
        index += 1

    return case_dict


def create_dir(i, sign, index, no_element, axis):

	"""Creation of folders."""

	os.mkdir(f'Coord-{index}-{no_element}-{axis}-{i}_{sign}')

def rename(filename, name):

	"""Rename files."""

	os.rename(filename, name)


def move_file(src_path, dst_path, filename):

	"""Move files into directories."""

	src_path = os.path.join(src_path, filename)
	dst_path = os.path.join(dst_path, filename)
	shutil.move(src_path, dst_path)


def iterate_over_files(pattern, name, input_file_name, code):

	"""Iterate over all input files and move them 
	into subdirectories. Rename the input file names.
	"""

	drct = os.getcwd()
	files = Path(drct).glob(pattern)

	disp = displacements(input_file_name, code)

	for file in files:
		f = file
		f = str(f)
		namefile = f.split("/")
		namefile = namefile[-1]
		number = re.sub('[^0-9]', '', namefile)
		number = int(number)
		sign = even_odd_numbers(number)		

		for key, values in disp.items():
			if number == values[0] or number == values[1]:
				
				create_dir(namefile, sign, key, values[2], values[3])
				created_dir = f'Coord-{key}-{values[2]}-{values[3]}-{namefile}_{sign}'

				dst_path = os.path.join(drct, created_dir)		

				move_file(drct, dst_path, namefile)

				os.chdir(created_dir)
				rename(namefile, name)	
				os.chdir('..')
