'''It checks if the single point calculations have been succesfully finished.'''

#AUTHOR: Ariadni Boziki

import os
import re


def single_successful_output(output, code):

    '''It checks if the cell and geometry optimization calculation has been successful (FHIaims).'''

    path = os.getcwd()

    if not os.path.exists(output):
        raise FileNotFoundError(f'The output file {output} does not exist')

    if os.path.exists(output):

        with open(output) as f:
            if (code == 'aims'):
                pattern = 'Have a nice day.'
                if not pattern in f.read():
                    print(f'THE CELL AND GEOMETRY OPTIMIZATION CALCULATION WAS NOT SUCCESSFUL \ PROBLEM DURING THE FHIAIMS CALCULATION\n')
                    flag_exit = True
                else:

                    print(f'THE CELL AND GEOMETRY OPTIMIZATION CALCULATION WAS SUCCESSFUL\n')
                    print(f'*************************************************************************************************************************************************************\n')
                    flag_exit = False

            elif (code == 'dftb+'):
                pattern = 'DFTB+ running times'
                if not pattern in f.read():
                    print(f'THE CELL AND GEOMETRY OPTIMIZATION CALCULATION WAS NOT SUCCESSFUL \ PROBLEM DURING THE DFTB+ CALCULATION\n')
                    flag_exit = True
                else:

                    print(f'THE CELL AND GEOMETRY OPTIMIZATION CALCULATION WAS SUCCESSFUL\n')
                    print(f'*************************************************************************************************************************************************************\n')
                    flag_exit = False

    return flag_exit

def successful_output(output, pattern): 

    '''It checks if the single point calculation
    related to frozen phonon approximation
    has been successful.'''

    path = os.getcwd()
    contents = filter(os.path.isdir, os.listdir(os.getcwd()))
    for i in contents:
        if 'Coord' in i:
            drct = i
        path_drct = os.path.join(path, drct)
        os.chdir(path_drct)

        with open(output) as f:
            if not pattern in f.read():
                print(f'{drct} has not successfully run')

        os.chdir('../')


def successful_output_dispersion_DFTBplus(output, pattern1, pattern2):

    '''It checks if the single point calculation
    related to frozen phonon approximation as well as 
    polarizability calculation in DFTB+ has been successful
    while dispersion is applied.'''

    path = os.getcwd()
    contents = filter(os.path.isdir, os.listdir(os.getcwd()))
    for i in contents:
        if 'Coord' in i:
            drct = i
        path_drct = os.path.join(path, drct)
        os.chdir(path_drct)

        try:
            with open(output) as f:
                if not pattern1 in f.read():
                    print(f'{drct} has not successfully run')
        except IOError:
            print(f'Output does not exist in {drct}')

        os.chdir('polarizability')
        try:
            with open(output) as ff:
                if not pattern2 in ff.read():
                    print(f'{drct} polarizability has not successfully run')
        except IOError:
            print(f'Output does not exist {drct} / polarizability')

        os.chdir('../../')
