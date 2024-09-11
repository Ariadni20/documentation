'''Preparation of inputs for cell and geometry optimization / single point calculations'''

#AUTHOR: Ariadni Boziki

import os
import sys
from theseus import InputsPreparation as inputs
from theseus import Constants as c
import hsd


class InputsGenerator:

    '''Generation of FHIaims and DFTB+ inputs.

    Attributes
    ----------

    code : str
        Code
    kpoints: str 
        K-points (Three values: they are converted to integer)
    functional: str
        Functional (FHIaims)
    eev: str
        Convergence criterion for the self-consistency cycle, based on the sum of eigenvalues (FHIaims)
    rho: str
        Convergence criterion for the self-consistency cycle, based on the charge density (FHIaims)
    etot: str
        Convergence criterion for the self-consistency cycle, based on the total energy (FHIaims)
    forces: str
        Convergence criterion for the self-consistency cycle, based on energy derivatives/forces (FHIaims)
    sc_iter_limit: str
        Maximum number of s.c.f. cycles before a calculation is considered and abandoned (FHIaims)
    species: str
        Species defaults settings; basis set, integration grids, accuracy of the Hartree potential (FHIaims)
    frequencies: bool
    geometry: str
        Maximum residual force component per atom, below which the geometry relaxation is considered converged (geometry relaxation) (FHIaims)
    energy: str
        Energy amount by which a relaxation step can move upwards and is still accepted (FHIaims)
    steps: str
        Maximum number of steps after which a structure optimization will be aborted (FHIaims)
    pol_grid: str
        Polarization grid
    max_force_component: str / converted to float
        Optimization is stopped, if the force component with the maximal absolute value goes below this value (DFTB+)
    max_steps: str
        Maximum number of steps after which the optimization should stop (DFTB+)
    SCC_tolerance: str
        Stopping criteria for the scc. Tolerance for the maximum difference in any charge between two scc cycles (DFTB+)
    max_SCC_iterations: str
        Maximal number of scc cycles to reach convergence (DFTB+)
    output_file: str
        Output file
    dispersion: bool
    path: str
        Path of the directory from which the code has been executed
    package_path: str
        Path of the directory where the code is located
    '''

    def __init__(
        self, 
        code,
        kpoints,
        functional=None,
        eev=None,
        rho=None,
        etot=None,
        forces=None,
        sc_iter_limit=None,
        species=None,
        frequencies=None,
        geometry=None,
        energy=None,
        steps=None,
        pol_grid=None,
        max_force_component=None,
        max_steps=None,
        SCC_tolerance=None,
        max_SCC_iterations=None,
        output_file=None,
        dispersion=False
    ):

        self.code = code
        self.kpoints = kpoints
        self.functional = functional
        self.eev = eev
        self.rho = rho
        self.etot = etot
        self.forces = forces
        self.sc_iter_limit = sc_iter_limit
        self.species = species
        self.frequencies = frequencies
        self.geometry = geometry
        self.energy = energy
        self.steps = steps
        self.pol_grid = pol_grid
        self.max_force_component = max_force_component
        self.max_steps = max_steps
        self.SCC_tolerance = SCC_tolerance
        self.max_SCC_iterations = max_SCC_iterations
        self.output_file = output_file
        self.dispersion = dispersion
        self.path = os.getcwd()
        self.species_path = os.path.abspath(__file__)
        self.package_path = os.path.dirname(self.species_path)


    def FHIaims_control_file(self):

        '''Generation of the control file of FHIaims
    
        number_atype: int, list
            The number of each atom within the crystal structure
        atype: str, list
            The atom type of each atom within the crystal structure
        unique_atype: str, list
            The atom types present in the crystal structure (appearing only once)
        input_parameters: dict
            It stores the parameters found in the control file and is utilized to display the input information in the output
        species_path_updated: str
            Directory path where the species defaults settings are located
        path_of_file_atype: str
            Path of the species defaults of an atom type
        source_data: str
            Data of the species defaults of an atom type
        '''

        number_atype, atype = inputs.read_the_atom_type('geometry.in', self.code)

        unique_atype = set()
        for string in atype:
            unique_atype.update(string)

        control_file = os.path.join(self.path, "control.in")

        with open(control_file, 'w') as fh:
            lines = [
                f'xc {self.functional}\n',
                f'k_grid {self.kpoints}\n',
                f'k_offset 0.5 0.5 0.5\n',
                f'relativistic atomic_zora scalar\n',
                f'\n',
                f'sc_accuracy_eev {self.eev}\n',
                f'sc_accuracy_rho {self.rho}\n',
                f'sc_accuracy_etot {self.etot}\n',
                f'sc_accuracy_forces {self.forces}\n',
                f'sc_iter_limit {self.sc_iter_limit}\n',
                f'\n'
            ]

            
            if self.dispersion:
                lines.extend([
                    f'many_body_dispersion_nl   beta=0.83\n'
                    ])


            if self.geometry is not None:
                lines.extend([
                    f'relax_geometry bfgs {self.geometry}\n',
                    f'energy_tolerance {self.energy}\n',
                    f'max_relaxation_steps {self.steps}\n',
                    f'relax_unit_cell full\n'
                    f'\n'
                ])

                input_parameters = {
                'FUNCTIONAL:': self.functional,
                'K-POINTS:': self.kpoints,
                'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON THE SUM OF EIGENVALUES:': self.eev,
                'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON THE CHARGE DENSITY:': self.rho,
                'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON THE TOTAL ENERGY:': self.etot,
                'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON ENERGY DERIVATIVES ("FORCES"):': self.forces,
                'MAXIMUM NUMBER OF S.C.F. CYCLES BEFORE A CALCULATION IS CONSIDERED AND ABANDONED:': self.sc_iter_limit,
                'SPECIES DEFAULTS SETTINGS (BASIS SET, INTEGRATION GRIDS, ACCURACY OF THE HARTREE POTENTIAL):': self.species,
                'MAXIMUM RESIDUAL FORCE COMPONENT PER ATOM (in eV/Å) BELOW WHICH THE GEOMETRY RELAXATION IS CONSIDERED CONVERGED (GEOMETRY RELAXATION):': self.geometry,
                'ENERGY AMOUNT BY WHICH A RELAXATION STEP CAN MOVE UPWARDS AND IS STILL ACCEPTED:': self.energy,
                'MAXIMUM NUMBER OF STEPS AFTER WHICH A STRUCTURE OPTIMIZATION WILL BE ABORTED:': self.steps
                }

                print(f'***************************************************')
                print(f'* Cell and geometry optimization input parameters *')
                print(f'***************************************************')
                for key, value in input_parameters.items():
                    print(f' {key} {value} ')
                print(f'\n')

            if self.frequencies:
                lines.extend([
                    f'compute_forces .true.\n',
                    f'final_forces_cleaned .true.\n'
                ])


                input_parameters = {
                    'FUNCTIONAL:': self.functional,
                    'K-POINTS:': self.kpoints,
                    'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON THE SUM OF EIGENVALUES:': self.eev,
                    'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON THE CHARGE DENSITY:': self.rho,
                    'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON THE TOTAL ENERGY:': self.etot,
                    'CONVERGENCE CRITERION FOR THE SELF-CONSISTENCY CYCLE, BASED ON ENERGY DERIVATIVES ("FORCES"):': self.forces,
                    'MAXIMUM NUMBER OF S.C.F. CYCLES BEFORE A CALCULATION IS CONSIDERED AND ABANDONED:': self.sc_iter_limit,
                    'SPECIES DEFAULTS SETTINGS (BASIS SET, INTEGRATION GRIDS, ACCURACY OF THE HARTREE POTENTIAL):': self.species,
                    'POLARIZATION GRID:': self.pol_grid
                }

                print(f'************************************************************************************************')
                print(f'* Single point calculation input parameters - Forces - Polarizability - Cartesian Polarization *')
                print(f'************************************************************************************************')
                for key, value in input_parameters.items():
                    print(f' {key} {value} ')
                print(f'\n')


            if self.pol_grid is not None:
                grid0 = self.pol_grid.split()[0]
                grid1 = self.pol_grid.split()[1]
                grid2 = self.pol_grid.split()[2]

                lines.extend([
                    f'DFPT dielectric\n',
                    f'KS_method serial\n',
                    f'output polarization 1 {grid0} 1 1\n',
                    f'output polarization 2 1 {grid1} 1\n',
                    f'output polarization 3 1 1 {grid2}\n',
                    '\n'
                    ])

            fh.writelines(lines)

        species_path_updated = os.path.join(self.package_path, 'species_defaults', self.species)

        for i in unique_atype:
            for filename in os.listdir(species_path_updated):
                species_type = filename.split('_')[1]
                if (i==species_type):
                    path_of_file_atype = os.path.join(species_path_updated, filename)
                    with open(path_of_file_atype, 'r') as source_file:
                        source_data = source_file.read()
                    with open(control_file, 'a') as target_file:
                        target_file.write(source_data)



    def _dict_to_hsd(self, data):

        '''It returns the necessary information in the correct format for DFTB+ to read the data specified in GenFormat/Geometry.
        '''

        hsd = ""
        for key, value in data.items():
            if isinstance(value, dict):
                if key == "Geometry":            
                    hsd += f"{key} {{"
                    hsd += self._dict_to_hsd(value)
            else:
                if key == "GenFormat":
                    hsd += f"{key} {{\n"
                    hsd += f" <<<'{value}'\n"
                    hsd += "}}\n"
        return hsd



    def DFTB_parameters_input_file(self):

        '''It generates the dftb_in.hsd input file of DFTB+

        species_path_updated: str
            Directory path where the .skf files are located
        data: dict
            It stores the parameters found in the dftb_in.hsd file
        '''

        species_path_updated = os.path.join(self.package_path, '3ob-3-1')

        number_atom_type, atom_type = inputs.read_the_atom_type('geo.gen', 'dftb+')

        atype_unique_set = set(atom_type)
        atype_unique_list = list(atype_unique_set)

        kpoint0 = self.kpoints.split()[0]
        kpoint1 = self.kpoints.split()[1]
        kpoint2 = self.kpoints.split()[2]

        value0 = 1 if int(kpoint0) % 2 == 0 else 0
        value1 = 1 if int(kpoint1) % 2 == 0 else 0
        value2 = 1 if int(kpoint2) % 2 == 0 else 0

        if self.max_force_component is not None:

            if self.dispersion:

                data = {
                    'Geometry': {
                        "GenFormat": "geo.gen"
                    },
                    'Driver': 
                        {'ConjugateGradient': 
                            {'MovedAtoms': '1:-1', 
                            'MaxForceComponent': float(self.max_force_component), 
                            'MaxSteps': self.max_steps, 
                            'OutputPrefix': f'"{self.output_file}"', 
                            'LatticeOpt': True}}, 
                    'Hamiltonian': 
                        {'DFTB': {'SCC': True, 
                            'SCCTolerance': self.SCC_tolerance, 
                            'MaxSCCIterations': self.max_SCC_iterations, 
                            'SlaterKosterFiles': {'Type2FileNames': 
                                {'Prefix': f'"{species_path_updated}/"', 'Separator': '"-"', 'Suffix': '".skf"'}}, 
                            'MaxAngularMomentum': {i: f'"{c.max_angular_mom[i]}"' for i in atype_unique_list},
                            'ThirdOrderFull': True, 
                            'HCorrection': {'Damping': {'Exponent': 4.05}}, 
                            'HubbardDerivs': {i: c.hubbard_derivates[i] for i in atype_unique_list}, 
                            'KPointsAndWeights': {'SuperCellFolding': [[kpoint0, 0, 0], [0, kpoint1, 0], [0, 0, kpoint2], [value0, value1, value2]]},
                            'Dispersion': {'MBD': {'KGrid': [1, 1, 1], 'Beta': 0.83}}}},
                    'Options': {'WriteResultsTag': True}, 
                    'Analysis': {'CalculateForces': True}} 

                input_parameters_dftb = {
                    'OPTIMIZATION IS STOPPED, IF THE FORCE COMPONENT WITH THE MAXIMAL ABSOLUTE VALUE GOES BELOW:': self.max_force_component,
                    'MAXIMUM NUMBER OF STEPS AFTER WHICH THE OPTIMIZATION SHOULD STOP:': self.max_steps,
                    'STOPPING CRITERIA FOR THE SCC. TOLERANCE FOR THE MAXIMUM DIFFERENCE IN ANY CHARGE BETWEEN TWO SCC CYCLES:': self.SCC_tolerance,
                    'MAXIMAL NUMBER OF SCC CYCLES TO REACH CONVERGENCE:': self.max_SCC_iterations,
                    'K-POINTS:': self.kpoints,
                    'DISPERSION:': 'MBD'
                        }

                print(f'***************************************************')
                print(f'* Cell and geometry optimization input parameters *')
                print(f'***************************************************')
                for key, value in input_parameters_dftb.items():
                    print(f' {key} {value} ')
                print(f'\n')


            else:

                data = {
                    'Geometry': {
                        "GenFormat": "geo.gen"
                    },
                    'Driver': 
                        {'ConjugateGradient': 
                            {'MovedAtoms': '1:-1', 
                            'MaxForceComponent': float(self.max_force_component), 
                            'MaxSteps': self.max_steps, 
                            'OutputPrefix': f'"{self.output_file}"', 
                            'LatticeOpt': True}}, 
                    'Hamiltonian': 
                        {'DFTB': {'SCC': True, 
                            'SCCTolerance': self.SCC_tolerance, 
                            'MaxSCCIterations': self.max_SCC_iterations, 
                            'SlaterKosterFiles': {'Type2FileNames': 
                                {'Prefix': f'"{species_path_updated}/"', 'Separator': '"-"', 'Suffix': '".skf"'}}, 
                            'MaxAngularMomentum': {i: f'"{c.max_angular_mom[i]}"' for i in atype_unique_list},
                            'ThirdOrderFull': True, 
                            'HCorrection': {'Damping': {'Exponent': 4.05}}, 
                            'HubbardDerivs': {i: c.hubbard_derivates[i] for i in atype_unique_list}, 
                            'KPointsAndWeights': {'SuperCellFolding': [[kpoint0, 0, 0], [0, kpoint1, 0], [0, 0, kpoint2], [value0, value1, value2]]}}}, 
                    'Options': {'WriteResultsTag': True}, 
                    'Analysis': {'CalculateForces': True}} 

                input_parameters_dftb = {
                    'OPTIMIZATION IS STOPPED, IF THE FORCE COMPONENT WITH THE MAXIMAL ABSOLUTE VALUE GOES BELOW:': self.max_force_component,
                    'MAXIMUM NUMBER OF STEPS AFTER WHICH THE OPTIMIZATION SHOULD STOP:': self.max_steps,
                    'STOPPING CRITERIA FOR THE SCC. TOLERANCE FOR THE MAXIMUM DIFFERENCE IN ANY CHARGE BETWEEN TWO SCC CYCLES:': self.SCC_tolerance,
                    'MAXIMAL NUMBER OF SCC CYCLES TO REACH CONVERGENCE:': self.max_SCC_iterations,
                    'K-POINTS:': self.kpoints
                        }

                print(f'***************************************************')
                print(f'* Cell and geometry optimization input parameters *')
                print(f'***************************************************')
                for key, value in input_parameters_dftb.items():
                    print(f' {key} {value} ')
                print(f'\n')

        else:

            if self.dispersion:

                data = {
                    'Geometry': {
                        "GenFormat": "geo.gen"
                    },
                    'Hamiltonian': 
                        {'DFTB': {'SCC': True, 
                            'SCCTolerance': self.SCC_tolerance, 
                            'MaxSCCIterations': self.max_SCC_iterations, 
                            'SlaterKosterFiles': {'Type2FileNames': 
                                {'Prefix': f'"{species_path_updated}/"', 'Separator': '"-"', 'Suffix': '".skf"'}}, 
                            'MaxAngularMomentum': {i: f'"{c.max_angular_mom[i]}"' for i in atype_unique_list},
                            'ThirdOrderFull': True, 
                            'HCorrection': {'Damping': {'Exponent': 4.05}}, 
                            'HubbardDerivs': {i: c.hubbard_derivates[i] for i in atype_unique_list}, 
                            'KPointsAndWeights': {'SuperCellFolding': [[kpoint0, 0, 0], [0, kpoint1, 0], [0, 0, kpoint2], [value0, value1, value2]]},
                            'Dispersion': {'MBD': {'KGrid': [1, 1, 1], 'Beta': 0.83}}}},
                    'Options': {'WriteResultsTag': True}, 
                    'Analysis': {'CalculateForces': True}}

                input_parameters_dftb = {
                    'STOPPING CRITERIA FOR THE SCC. TOLERANCE FOR THE MAXIMUM DIFFERENCE IN ANY CHARGE BETWEEN TWO SCC CYCLES:': self.SCC_tolerance,
                    'MAXIMAL NUMBER OF SCC CYCLES TO REACH CONVERGENCE:': self.max_SCC_iterations,
                    'K-POINTS:': self.kpoints,
                    'DISPERSION:': 'MBD'
                        }

                print(f'****************************')
                print(f'* Single point calculation *')
                print(f'****************************')

                for key, value in input_parameters_dftb.items():
                    print(f' {key} {value} ')
                print(f'\n')

            else:

                data = {
                    'Geometry': {
                        "GenFormat": "geo.gen"
                    },
                    'Hamiltonian': 
                        {'DFTB': {'SCC': True, 
                            'SCCTolerance': self.SCC_tolerance, 
                            'MaxSCCIterations': self.max_SCC_iterations, 
                            'SlaterKosterFiles': {'Type2FileNames': 
                                {'Prefix': f'"{species_path_updated}/"', 'Separator': '"-"', 'Suffix': '".skf"'}}, 
                            'MaxAngularMomentum': {i: f'"{c.max_angular_mom[i]}"' for i in atype_unique_list},
                            'ThirdOrderFull': True, 
                            'HCorrection': {'Damping': {'Exponent': 4.05}}, 
                            'HubbardDerivs': {i: c.hubbard_derivates[i] for i in atype_unique_list}, 
                            'KPointsAndWeights': {'SuperCellFolding': [[kpoint0, 0, 0], [0, kpoint1, 0], [0, 0, kpoint2], [value0, value1, value2]]}}}, 
                    'Options': {'WriteResultsTag': True}, 
                    'Analysis': {'CalculateForces': True, 'Polarisability': {}, 'DegeneracyTolerance' : 1024}}

                input_parameters_dftb = {
                    'STOPPING CRITERIA FOR THE SCC. TOLERANCE FOR THE MAXIMUM DIFFERENCE IN ANY CHARGE BETWEEN TWO SCC CYCLES:': self.SCC_tolerance,
                    'MAXIMAL NUMBER OF SCC CYCLES TO REACH CONVERGENCE:': self.max_SCC_iterations,
                    'K-POINTS:': self.kpoints
                        }

                print(f'****************************')
                print(f'* Single point calculation *')
                print(f'****************************')

                for key, value in input_parameters_dftb.items():
                    print(f' {key} {value} ')
                print(f'\n')


        hsd_content = self._dict_to_hsd(data)
        data_without_genformat = data.copy()
        data_without_genformat["Geometry"].pop("GenFormat")
        data_without_genformat.pop("Geometry")

        hsd.dump(data_without_genformat, "output1.hsd")
        with open("dftb_in.hsd", "w") as hsd_file:
            hsd_file.write(hsd_content)

        with open("output1.hsd", "r") as hsd_file1:
            content_to_append = hsd_file1.read()

        with open("dftb_in.hsd", "a") as hsd_file:
            hsd_file.write(content_to_append)

        os.remove("output1.hsd")
