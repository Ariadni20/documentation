'''Submission of FHI-aims and PHONOPY calculations'''

#AUTHOR: Ariadni Boziki

import os
import re
import subprocess
import shutil
import concurrent.futures
import numpy as np
from theseus import InputsPreparation as inputs
from theseus import FiniteDisplacements_phonopy as finitedisps
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from phonopy import Phonopy
from phonopy.interface.calculator import read_crystal_structure, write_crystal_structure, write_supercells_with_displacements, get_default_displacement_distance
from phonopy.file_IO import write_FORCE_SETS, write_FORCE_CONSTANTS
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections



class PhonopyCalculator:
    
    def __init__(self, code, cell_dims, output_file_of_SCFSs):
        
        self.code = code
        self.cell_dims = cell_dims
        self.output_file_of_SCFSs = output_file_of_SCFSs

    
    def _read_crystal_structure(self):

        if (self.code == 'aims'):
            return read_crystal_structure('geometry.in', interface_mode='aims')
        elif (self.code == 'dftb+'):
            return read_crystal_structure('geo.gen', interface_mode='dftbp')


    def _setup_phonopy(self):

        unitcell, optional_structure_info = self._read_crystal_structure()
        cell_dims_int = [int(num) for num in self.cell_dims.split()]
        supercell_matrix = [[cell_dims_int[0], 0, 0], [0, cell_dims_int[1], 0], [0, 0, cell_dims_int[2]]]
        self.phonon = Phonopy(unitcell,supercell_matrix)
        default_displacement_for_code = get_default_displacement_distance(interface_mode=self.code)
        self.phonon.generate_displacements(distance=default_displacement_for_code)
        self.disps = self.phonon.displacements


    def supercell_disp_PHONOPY(self):

        self._setup_phonopy()
        supercells = self.phonon.supercells_with_displacements
        gen_supercell = self.phonon.supercell
        if (self.code == 'aims'):
            write_supercells_with_displacements('aims', gen_supercell, supercells)
        elif (self.code == 'dftb+'):
            write_supercells_with_displacements('dftbp', gen_supercell, supercells)


    def _number_of_atoms(self):

        if (self.code == 'aims'):
            geometry_processor = inputs.GeometryProcessor('geometry.in.supercell', self.code)
        elif (self.code == 'dftb+'):
            geometry_processor = inputs.GeometryProcessor('geo.genS', self.code)
        
        self.no_of_atoms = geometry_processor.number_of_atoms()
            

    def _read_forces_from_output(self, output_path):

        '''Returns the forces from the FHIaims output and results.tag DFTB+ file.'''

        self._number_of_atoms()

        if (self.code == 'aims'):
            forces_identifier = 'Total atomic forces (unitary forces cleaned)'
        elif (self.code == 'dftb+'):
            forces_identifier = 'forces'

        f=open(output_path, 'r')
        found_unique_line = False

        forces_lines = []
        for line in f:
            if forces_identifier in line:
                found_unique_line = True
                continue

            if found_unique_line:
                forces_lines.append(line.strip())
                if len(forces_lines) >= self.no_of_atoms:
                    break

        self.forces = []
        pattern = r"[-+]?\d+\.\d+E[-+]?\d+"

        for item in forces_lines:
            numbers = re.findall(pattern, item)
            if len(numbers) >= 3:
                last_three_columns = [float(num) for num in numbers[-3:]]
                self.forces.append(last_three_columns)

        return self.forces


    def disp_forces_dataset_dyn_matrix(self):

        '''Generates the FORCE_SETS file of PHONOPY'''

        self._number_of_atoms()
        self.dataset = {'natom': self.no_of_atoms,
                'first_atoms': []}

        path = os.getcwd()
     
        contents = filter(os.path.isdir, os.listdir(os.getcwd()))
        conts = []
        for ii in contents:
            if 'Coord' in ii:
                conts.append(ii)
        conts.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

        self._setup_phonopy()

        for j, jj in zip(self.disps,conts):

            if (self.code == 'aims'):
                output_path = os.path.join(path, jj, self.output_file_of_SCFSs)
            elif (self.code == 'dftb+'):
                output_path = os.path.join(path, jj, 'results.tag')

            self.forces = self._read_forces_from_output(output_path)

            entry = {
                    'number': np.array([j[0]]),
                    'displacement': np.array(j[1:]),
                    'forces': np.array(self.forces)
            }
            self.dataset['first_atoms'].append(entry)

        self.phonon.dataset = self.dataset

        self.phonon.produce_force_constants()

        q = [0.0,0.0,0.0]
        self.dynamical_matrix = self.phonon.get_dynamical_matrix_at_q(q)

        return self.dataset, self.dynamical_matrix


    def plot_band_structure(self):


        path = [[[0, 0, 0], [0, 0.5, 0], [0, 0.5, 0.5]],
                [[0, 0, 0.5], [0, 0, 0], [-0.5, 0, 0.5], [-0.5, 0.5, 0.5]],
                [[0, 0.5, 0], [-0.5, 0.5, 0], [-0.5, 0, 0], [0, 0, 0]]]
        labels = ["$\\Gamma$", "Z", "D", "B", "$\\Gamma$", "A", "E", "Z", "$C_{2}$", "$Y_{2}$", "$\\Gamma$"]
        qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
        self.phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
        self.phonon.plot_band_structure().savefig("phonon_band_structure.pdf")
        self.phonon.write_yaml_band_structure()
        band_dict = self.phonon.get_band_structure_dict()



class Calculator:


    def __init__(self, code, output_file, commands=None):
    
        self.code = code
        self.commands = commands
        self.output_file = output_file
        contents = filter(os.path.isdir, os.listdir(os.getcwd()))
        conts = []
        for ii in contents:
            if 'Coord' in ii:
                conts.append(ii)
        conts.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
        self.cont = conts


    def submit_job(self):

        '''Submission of FHIaims calculations'''

        subprocess.run(self.commands, shell=True, executable='/bin/bash')


    def frozen_phonon_approximation_drct(self):

        '''Creation of the directory where the files associated with 
            frozen phonon approximation will be saved.'''

        if (self.code == 'aims'):
            geometry = 'geometry.in'
            new_geometry = 'geometry.in.next_step'
        elif (self.code == 'dftb+'):
            geometry = 'geo.gen'
            new_geometry = f'{self.output_file}.gen'

        path = os.getcwd()
        os.mkdir(f'vibrations')

        output_path = os.path.join(path, 'vibrations')

        geometry_path = os.path.join(path, new_geometry)

        if os.path.exists(geometry_path):
            shutil.copy(new_geometry, output_path)  
            
        else:
            shutil.copy(geometry, output_path)  
      
        os.chdir(output_path)

        geometry_path = os.path.join(output_path, new_geometry)

        if os.path.exists(geometry_path):
            os.rename(new_geometry, geometry) 


    def run_command(self, command):
        try:
            # Run the command and capture its output
            result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
            # Process the result as needed
            return result.stdout
        except subprocess.CalledProcessError as e:
            # Handle command failure
            return f"Error running {command}: {e}"


    def submit_jobs_in_parallel(self):

        # Define a list of commands to run in parallel


        command_statement = []

        for i in self.cont:
    
            command_tmp = f'cd {i}; {self.commands}'
            command_statement.append(command_tmp)

        number_of_nodes = 1
        number_of_cores = 128

        no_of_folders = len(command_statement)

        # Create a ThreadPoolExecutor with a specified number of worker threads
        num_threads = int((number_of_nodes*number_of_cores)/no_of_folders)  # Adjust the number of threads as needed
        print(num_threads)
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            # Submit the commands for execution in parallel
            futures = [executor.submit(self.run_command, command) for command in command_statement]
            # Wait for all commands to complete and get their results
            results = [future.result() for future in futures]
        # Process the results
        for i, result in enumerate(results):
            print(f"Result for command {i + 1}:\n{result}")



    def run_parallel_via_slurm(self):


        command_statement = []

        for i in self.cont:
    
            command_tmp = f'cd {i}; {self.commands}'
            command_statement.append(command_tmp)

        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(self.run_command, command_statement)
