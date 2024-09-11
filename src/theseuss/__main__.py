'''Main module of the package'''

#AUTHOR: Ariadni Boziki


import os
import subprocess
from theseus import InputOutputFiles as inputfiles
from theseus import FiniteDisplacements_phonopy as disp
from theseus import TwoPointCentralDifference as cendiff
from theseus import EigenvectorsFrequenciesPHONOPY as eigvec
from theseus import IRRamanSpectra as spectra
from theseus import InputsPreparation as inputs
from theseus import SubmitCalculations as submit
from theseus import MapAtoms
from theseus import PlotSpectra as plots
from theseus import Coordinates as coords
from theseus import CheckSuccessOutput as check
import numpy as np
import argparse
import re



def main(cell_geometry_optimization, functional, kpoints, eev, rho, etot, forces, sc_iter_limit, species, geometry, energy, steps, pol_grid, supercell, code, output_file, dispersion, spectra_calculation, files_preparation, phonons, cell_dims, submission_cell, plot_bands, commands, max_force_component, max_steps, SCC_tolerance, max_SCC_iterations):


    message = 'Input'
    border = '*' * (len(message) + 4)
    print(f'{border}\n* {message} *\n{border}')
    print(f'CELL AND GEOMETRY OPTIMIZATION: {cell_geometry_optimization}')
    print(f'SUPERCELL: {supercell}')
    print(f'DISPERSION: {dispersion}')
    print(f'SUBMISSION OF PHONOPY JOB TO GENERATE DISPLACEMENTS: {submission_cell}')
    print(f'PREPARATION OF FILES FOR FROZEN-PHONON APPROXIMATION: {files_preparation}')
    print(f'PHONONS CALCULATION: {phonons}')
    print(f'CALCULATION OF SPECTRA: {spectra_calculation}')
    print(f'DIMENSIONS OF SUPERCELL: {cell_dims}')
    print(f'CODE: {code}')
    print(f'OUTPUT FILE: {output_file}\n')


    calculator = submit.Calculator(code, output_file, commands)
    initial_space_group = MapAtoms.get_international_space_group_number(code)
    print(f'INTERNATIONAL SPACE GROUP NUMBER OF THE EXPERIMENTAL STRUCTURE: {initial_space_group}\n')


    if cell_geometry_optimization:

        generator = inputfiles.InputsGenerator(code, kpoints, functional, eev, rho, etot, forces, sc_iter_limit, species,
            False, geometry, energy, steps, None, max_force_component, max_steps, SCC_tolerance, max_SCC_iterations, 
            output_file, dispersion)

        if (code == 'aims'):

            if os.path.isfile('control.in'):
                print(f'THE CONTROL FILE (input of FHI-aims) FOR CELL AND GEOMETRY OPTIMIZATION ALREADY EXISTS\n')
                pass
            else:
                try:
                    generator.FHIaims_control_file()
                    print(f'THE CONTROL FILE (input of FHI-aims) FOR CELL AND GEOMETRY OPTIMIZATION HAS BEEN GENERATED\n')
                except Exception as e:
                    print(f'THE CONTROL FILE GENERATION WAS NOT SUCCESSFUL / AN ERROR WAS OCCURED: {e}')

        elif (code == 'dftb+'):

            if os.path.isfile('dftb_in.hsd'):
                print(f'THE dftb_in.hsd FILE (input of DFTB+) FOR CELL AND GEOMETRY OPTIMIZATION ALREADY EXISTS\n')
                pass
            else:
                try:
                    generator.DFTB_parameters_input_file()
                    print(f'THE dftb_in.hsd FILE (input of DFTB+) FOR CELL AND GEOMETRY OPTIMIZATION HAS BEEN GENERATED\n')
                except Exception as e:
                    print(f'THE dftb_in.hsd FILE GENERATION WAS NOT SUCCESSFUL / AN ERROR WAS OCCURED: {e}')

        calculator.submit_job()

        flag_exit = check.single_successful_output(output_file, code)
        if flag_exit:
            exit()
        else:
            calculator.frozen_phonon_approximation_drct()
            space_group_of_optimized_str = MapAtoms.get_international_space_group_number(code)
            print(f'INTERNATIONAL SPACE GROUP NUMBER OF THE OPTIMIZED STRUCTURE: {space_group_of_optimized_str}\n')


    if submission_cell:    

        try:
            phonopy_calculator = submit.PhonopyCalculator(code, cell_dims, output_file)
            phonopy_calculator.supercell_disp_PHONOPY()

            print(f'PHONOPY GENERATED THE STRUCTURES WITH THE DISPLACEMENTS\n')
            print('*' * 150) 

        except Exception as e:
            print(f'GENERATION OF THE STRUCTURES WITH THE DISPLACEMENTS BY PHONOPY / AN ERROR WAS OCCURED: {e}')
            print('*' * 150) 


    if files_preparation:

        generator = inputfiles.InputsGenerator(code, kpoints, functional, eev, rho, etot, forces, sc_iter_limit, species, True, None, None, None, pol_grid, None, None, SCC_tolerance, max_SCC_iterations, output_file, dispersion)

        if (code == 'aims'):

            try:

                pattern = 'geometry.in-*'
                name = 'geometry.in'
                input_file_name = 'geometry.in-'
                disp.iterate_over_files(pattern, name, input_file_name, code)

                try:
                    generator.FHIaims_control_file()
                    print(f'THE CONTROL FILE (input of FHI-aims) FOR SINGLE POINT CALCULATION HAS BEEN GENERATED\n')
                except Exception as e:
                    print(f'THE CONTROL FILE (input of FHI-aims) FOR SINGLE POINT CALCULATION HAS NOT BEEN GENERATED / AN ERROR WAS OCCURED: {e}')

                inputs.copy_files_in_dir('Coord-*', 'control.in')

                print(f'THE DIRECTORIES THAT CONTAIN THE DISPLACED STRUCTURES HAVE BEEN GENERATED\n')
                print('*' * 150) 

            except Exception as e:

                print(f'THE DIRECTORIES THAT CONTAIN THE DISPLACED STRUCTURES HAVE NOT BEEN GENERATED / AN ERROR WAS OCCURED: {e}')


        if (code == 'dftb+'):

            try:

                pattern = 'geo.genS-*'
                name = 'geo.gen'
                input_file_name = 'geo.genS-'
                disp.iterate_over_files(pattern, name, input_file_name, code)

                try: 
                    generator.DFTB_parameters_input_file()
                    print(f'THE dftb_in.hsd (input of DFTB+) FOR SINGLE POINT CALCULATION HAS BEEN GENERATED\n')
                except Exception as e:
                    print(f'THE dftb_in.hsd FILE (input of DFTB+) FOR SINGLE POINT CALCULATION HAS NOT BEEN GENERATED / AN ERROR WAS OCCURED: {e}')

                inputs.copy_files_in_dir('Coord-*', 'dftb_in.hsd')

                if dispersion:

                    generator = inputfiles.InputsGenerator(code, kpoints, functional, eev, rho, etot, forces, 
                        sc_iter_limit, species, True, None, None, None, pol_grid, None, None, 
                        SCC_tolerance, max_SCC_iterations, output_file, None)
                    path = os.getcwd()
                    contents = filter(os.path.isdir, os.listdir(os.getcwd()))
                    for i in contents:
                        if 'Coord' in i:
                            drct = i
                        path_drct = os.path.join(path, drct)
                        os.chdir(path_drct)

                        os.mkdir(f'polarizability')	
                        inputs.copy_files_in_dir('polarizability', 'geo.gen')
                        os.chdir('polarizability')
                        generator.DFTB_parameters_input_file()
                        os.chdir('../../')	
    
                    print(f'POLARIZABILITY DIRECTORIES HAVE BEEN GENERATED WITHIN THE DIRECTORIES THAT INCLUDE THE DISPLACED STRUCTURES SINCE A DISPERSION METHOD IS USED\n')
                print(f'THE DIRECTORIES THAT CONTAIN THE DISPLACED STRUCTURES HAVE BEEN GENERATED\n')
                print('*' * 150)

            except Exception as e:

                print(f'THE DIRECTORIES THAT CONTAIN THE DISPLACED STRUCTURES HAVE NOT BEEN GENERATED / AN ERROR WAS OCCURED: {e}')

        FHIaims_calculator = submit.Calculator(code, output_file, commands)
#        FHIaims_calculator.submit_jobs_in_parallel()
        FHIaims_calculator.run_parallel_via_slurm()

    if phonons:

        try:

            phonopy_calculator = submit.PhonopyCalculator(code, cell_dims, output_file)
            dataset, dynamical_matrix = phonopy_calculator.disp_forces_dataset_dyn_matrix() 
            
            print(f'THE FORCE CONSTANTS HAVE BEEN CALCULATED USING PHONOPY\n')
            print(f'THE DYNAMICAL MATRIX HAS BEEN CALCULATED USING PHONOPY\n')
            print('*' * 150)

        except Exception as e:
            print(f'FORCE CONSTANTS/DYNAMICAL MATRIX / AN ERROR WAS OCCURED: {e}')
            print('*' * 150)


    if plot_bands:


        phonopy_calculator.plot_band_structure()


    if spectra_calculation:

        if (code == 'aims'):
            check.successful_output(output_file, 'Have a nice day.')
        if (code == 'dftb+') and not dispersion:
            check.successful_output(output_file, 'DFTB+ running times')
        if (code == 'dftb+') and dispersion:
            check.successful_output_dispersion_DFTBplus(output_file, 'DFTB+ running times', 'DFTB+ running times')

        phonopy_calculator = submit.PhonopyCalculator(code, cell_dims, output_file)
        dataset, dynamical_matrix = phonopy_calculator.disp_forces_dataset_dyn_matrix() 

        dynamical_matrix = np.real(dynamical_matrix)

        eigvecs, eigvals, freq = eigvec.read_eig_vec_phonopy(code, dynamical_matrix)
        print(f'eigvecs {eigvecs}')
        print(f'FREQUENCIES (1/cm)')
        print(f'{freq}')

        path = os.getcwd()
        contents = filter(os.path.isdir, os.listdir(os.getcwd()))
        for i in contents:
            if 'Coord' in i:
                drct = i
        path_drct = os.path.join(path, drct)
        os.chdir(path_drct)

        volume = cendiff.UnitCellVolume(output_file)
        print(f'VOLUME: {volume} A\N{SUPERSCRIPT THREE}')

        os.chdir('../')	

        sign_atom_coord = cendiff.subgroups_finite_disp_drct()
        pol, cartesian_pol, element_axis_coord = cendiff.pol_fhiaims(sign_atom_coord, 0.01, code, dispersion, output_file)


        if supercell:
            if (code == 'aims'):
                coords.cartesian_to_fractional_FHIaims(supercell)
            if (code == 'dftb+'):
                coords.cartesian_to_fractional_DFTBplus(supercell)
        else:
            if (code == 'aims'):
                with open('geometry.in','r') as file:
                    file_content = file.read()
                    if 'atom_frac' in file_content:
                        pass
                    else:
                        coords.cartesian_to_fractional_FHIaims(supercell)
            if (code == 'dftb+'):
                with open('geo.gen','r') as file:
                    file_content = file.read()
                    if 'F' in file_content:
                        pass
                    else:
                        coords.cartesian_to_fractional_DFTBplus(supercell)



        pol_extended, cartesian_pol_extended = MapAtoms.merge_arrays(pol, cartesian_pol, element_axis_coord)

        pol_ext, cartesian_pol_ext = MapAtoms.map_properties_of_equivalent_atoms(pol_extended, cartesian_pol_extended, code)

        pol, cartesian_pol = MapAtoms.remove_columns_from_array(pol_ext, cartesian_pol_ext)

###########################################################################
#   These commands are added if the cell has been multiplied.
        if supercell:
            if (code == 'aims'):
                geometry_processor = inputs.GeometryProcessor('geometry_backup.in', code) 
            elif (code == 'dftb+'):
                geometry_processor = inputs.GeometryProcessor('geo_backup.gen', code)
            
            number_of_atoms_unit_cell = geometry_processor.number_of_atoms()
            print(f'NUMBER OF ATOMS OF UNIT CELL: {number_of_atoms_unit_cell}')

            number_of_lines_unit_cell = 3*number_of_atoms_unit_cell

            if (code == 'aims'):
                geometry_processor = inputs.GeometryProcessor('geometry.in.supercell', code)
            elif (code == 'dftb+'):
                geometry_processor = inputs.GeometryProcessor('geo.genS', code)

            number_of_atoms = geometry_processor.number_of_atoms()
            print(f'NUMBER OF ATOMS OF SUPERCELL: {number_of_atoms}')

            number_of_lines = 3*number_of_atoms

            pol = pol[:number_of_lines_unit_cell]
            cartesian_pol = cartesian_pol[:number_of_lines_unit_cell]

###########################################################################

        print(f'THE SPECTRA HAVE BEEN CALCULATED\n')

        IRintensity = spectra.IR_intensity(volume, cartesian_pol, eigvecs)
        print(f'IR INTENSITY (D\N{SUPERSCRIPT TWO}\A\N{SUPERSCRIPT TWO} amu)')
        print(f'{IRintensity}')

        Ramanactivity = spectra.Raman_activity(pol, eigvecs, code)
        print(f'RAMAN ACTIVITY (A\N{SUPERSCRIPT FOUR}/amu)')
        print(f'{Ramanactivity}')

        np.savetxt("IRintensity.txt", IRintensity)
        np.savetxt("Ramanactivity.txt", Ramanactivity)
        np.savetxt("Frequency.txt", freq)

        plots.plot_spectrum_IR(freq, IRintensity, code)
        plots.plot_spectrum_Raman(freq, Ramanactivity, code)

        print('*' * 150) 


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', dest='input_file_name', type=str)
    args = parser.parse_args()
    input_file = args.input_file_name

    keywords = {} 

    with open(input_file, 'r') as file1:
        lines = file1.readlines()

        patt = r'(\S+)\s+(.+)'

        for line in lines:
            match = re.match(patt, line)
            if match:
                key = match.group(1)
                value = match.group(2)


                if key in ('cell_geometry_optimization', 'supercell', 'dispersion', 'spectra_calculation', 'files_preparation', 'phonons', 'submission', 'plot_bands'):
                    value = eval(value)
                if key in ('functional', 'eev', 'rho', 'etot', 'forces', 'sc_iter_limit', 'species', 'energy', 'geometry', 'steps', 'code', 'output_file', 'commands', 'max_force_component', 'max_steps', 'SCC_tolerance', 'max_SCC_iterations'):
                    value = value
                if key in ('kpoints', 'dimensions', 'pol_grid'):
                    value = ' '.join(value.split()[0:])

            keywords[key] = value


    cell_geometry_optimization = keywords.get('cell_geometry_optimization')
    supercell = keywords.get('supercell')
    dispersion = keywords.get('dispersion')
    spectra_calculation = keywords.get('spectra_calculation')
    files_preparation = keywords.get('files_preparation')
    phonons = keywords.get('phonons')
    submission_cell = keywords.get('submission')
    plot_bands = keywords.get('plot_bands')
    functional = keywords.get('functional')
    eev = keywords.get('eev')
    rho = keywords.get('rho')
    etot = keywords.get('etot')
    forces = keywords.get('forces')
    sc_iter_limit = keywords.get('sc_iter_limit')
    species = keywords.get('species')
    energy = keywords.get('energy')
    geometry = keywords.get('geometry')
    steps = keywords.get('steps')
    code = keywords.get('code')
    output_file = keywords.get('output_file')
    commands = keywords.get('commands')
    kpoints = keywords.get('kpoints')
    cell_dims = keywords.get('dimensions')
    pol_grid = keywords.get('pol_grid')
    max_force_component = keywords.get('max_force_component')
    max_steps = keywords.get('max_steps')
    SCC_tolerance = keywords.get('SCC_tolerance')
    max_SCC_iterations = keywords.get('max_SCC_iterations')


    main(cell_geometry_optimization, functional, kpoints, eev, rho, etot, forces, sc_iter_limit, species, geometry, energy, steps, pol_grid, supercell, code, output_file, dispersion, spectra_calculation, files_preparation, phonons, cell_dims, submission_cell, plot_bands, commands, max_force_component, max_steps, SCC_tolerance, max_SCC_iterations)

if __name__=='__main__':
    run()
