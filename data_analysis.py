# Module containing function which analyses parsed data from Gaussian/ORCA output files.

# Import modules


from math import exp
from rdkit.Chem import MolFromXYZBlock, rdMolTransforms  # Version 2022.9.4
import numpy as np  # Version 1.24.1
from scipy.optimize import linear_sum_assignment  # Version 1.14.1
from scipy.spatial import distance as scipy_distance


def analyse(parsed_data, settings):
    """Analyses parsed data from output files. Finds and excludes erroneous conformers and arranges data in a convenient
    format for writing to new files (e.g. supplementary information in a .docx file, where conformer data is arranged
    in order of increasing conformer energy). Also calculates Boltzmann-weighted averaged spectroscopic data."""

    # Define input data

    energies = parsed_data[0]
    element_list = parsed_data[1]
    x_cartesian_coords_list = parsed_data[2]
    y_cartesian_coords_list = parsed_data[3]
    z_cartesian_coords_list = parsed_data[4]
    wavelength_list = parsed_data[5]
    rotatory_strength_list = parsed_data[6]
    shielding_tensors = parsed_data[7]
    frequencies_list = parsed_data[8]
    sp_functional_and_basis_set = parsed_data[10]
    sp_solvent = parsed_data[11]
    sp_dispersion = parsed_data[12]
    opt_freq_functional_and_basis_set = parsed_data[13]
    opt_freq_solvent = parsed_data[14]
    opt_freq_dispersion = parsed_data[15]
    nmr_functional_and_basis_set = parsed_data[16]
    nmr_solvent = parsed_data[17]
    nmr_dispersion = parsed_data[18]
    tddft_functional_and_basis_set = parsed_data[19]
    tddft_solvent = parsed_data[20]
    tddft_dispersion = parsed_data[21]
    results_directory = parsed_data[23]
    oscillator_strength_list = parsed_data[24]
    list_of_conformer_suffixes = parsed_data[25]
    frequency_rotatory_strengths = parsed_data[29]
    optrot_wavelengths = parsed_data[30]
    optrot_strengths = parsed_data[31]
    frequency_dipole_strengths = parsed_data[32]
    or_functional_and_basis_set = parsed_data[33]
    or_solvent = parsed_data[34]
    or_dispersion = parsed_data[35]
    ir_intensities_list = parsed_data[36]
    chk_conf_suffixes = parsed_data[37]

    energy_threshold = float(settings["Energy cutoff (kcal/mol)"])
    mad_threshold = float(settings["MAD cutoff (A)"])
    c_slope = settings["C slope"]
    c_intercept = settings["C intercept"]
    h_slope = settings["H slope"]
    h_intercept = settings["H intercept"]
    temperature = float(settings["Temperature (K)"])

    # Define new lists

    ordered_energies = []
    ordered_relative_energies = []
    ordered_boltz_weights = []
    ordered_population_contributions = []
    ordered_element_list = []
    ordered_x_cartesian_coords_list = []
    ordered_y_cartesian_coords_list = []
    ordered_z_cartesian_coords_list = []
    ordered_frequencies_list = []
    ordered_ir_intensities_list = []
    ordered_frequency_rotatory_strengths_list = []
    ordered_frequency_dipole_strengths_list = []
    ordered_optrot_list = []
    ordered_optrot_wavelengths_list = []
    ordered_wavelength_list = []
    ordered_rotatory_strengths = []
    ordered_oscillator_strengths = []
    ordered_shielding_tensors = []
    duplicate_conformers = []
    dup_conf_numbers = []
    duplicate_conf_details = []
    relative_energies = []
    population_contributions = []
    boltz_weights = []
    boltz_c_tensors = []
    boltz_h_tensors = []
    boltz_shielding_tensors = []
    shielding_tensor_components = []
    frequency_components = []
    ir_intensity_components = []
    boltz_frequencies = []
    boltz_ir_intensities = []
    wavelength_components = []
    boltz_wavelengths = []
    rotatory_strength_components = []
    boltz_rotatory_strengths = []
    oscillator_strength_components = []
    boltz_oscillator_strengths = []
    frequency_rotatory_strength_components = []
    boltz_frequency_rotatory_strengths = []
    frequency_dipole_strength_components = []
    boltz_frequency_dipole_strengths = []
    optrot_strength_components = []
    boltz_optrot_strengths = []

    # Define variables

    hartree_to_kJmol = 2625.4996394799
    hartree_to_kcalmol = 627.5094740631
    boltz_constant_kJmol = 0.008314462618
    boltz_constant_kcalmol = 0.001987204259
    data_analysis_error_check = "No error detected"

    # Set appropriate Boltzmann constant

    if settings["Relative energy unit"] == "kJ/mol":
        boltz_constant = boltz_constant_kJmol
    elif settings["Relative energy unit"] == "kcal/mol":
        boltz_constant = boltz_constant_kcalmol
    # Remove conformers with imaginary frequency/ies, if present

    removals = 0
    for index in range(len(frequencies_list)):
        index = index - removals
        if min(frequencies_list[index]) < 0:
            # Remove non-converged conformer from data

            del energies[index]
            del element_list[index]
            del frequencies_list[index]
            del ir_intensities_list[index]
            del x_cartesian_coords_list[index]
            del y_cartesian_coords_list[index]
            del z_cartesian_coords_list[index]
            if list_of_conformer_suffixes:
                del list_of_conformer_suffixes[index]
            if chk_conf_suffixes:
                del chk_conf_suffixes[index]
            if wavelength_list:
                del wavelength_list[index]
                del rotatory_strength_list[index]
                del oscillator_strength_list[index]
            if shielding_tensors:
                del shielding_tensors[index]
            if frequency_rotatory_strengths:
                del frequency_rotatory_strengths[index]
                del frequency_dipole_strengths[index]
            if optrot_wavelengths:
                del optrot_strengths[index]
                del optrot_wavelengths[index]
            removals += 1  # This is because indices of remaining conformers are now reduced by 1, by removing the
            # last conformer
    number_imag_freq_confs_removed = removals

    # Check confs have same number of frequencies and IR intensities after removing confs with imaginary frequencies

    for conf in frequencies_list:
        a = len(frequencies_list[0])
        if len(conf) != a:
            data_analysis_error_check = "Error detected"
            error_message = "Inconsistent number of vibrational frequencies found\nfor conformers in file(s)."
            return (
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                data_analysis_error_check,
                error_message,
            )
    for conf in ir_intensities_list:
        a = len(ir_intensities_list[0])
        if len(conf) != a:
            data_analysis_error_check = "Error detected"
            error_message = "Inconsistent number of IR intensities found for conformers in file(s)."
            return (
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                "",
                data_analysis_error_check,
                error_message,
            )
    # Define function for building XYZ file format strings for a given conformer number, then creating an RDKit molecule

    def RDKitMoleculeMaker(conformer_number):
        """Function for building XYZ file format strings for a given conformer number,
        then creating an RDKit molecule."""
        xyz_block = ""
        for i in range(-2, len(element_list[conformer_number])):
            if i == -2:
                xyz_block += str(len(element_list[conformer_number]))
            elif i == -1:
                xyz_block += "\nMolecule"
            elif i > -1:
                if float(x_cartesian_coords_list[conformer_number][i]) >= 0:
                    x = "      " + str(x_cartesian_coords_list[conformer_number][i])
                else:
                    x = "     " + str(x_cartesian_coords_list[conformer_number][i])
                if float(y_cartesian_coords_list[conformer_number][i]) >= 0:
                    y = "      " + str(y_cartesian_coords_list[conformer_number][i])
                else:
                    y = "     " + str(y_cartesian_coords_list[conformer_number][i])
                if float(z_cartesian_coords_list[conformer_number][i]) >= 0:
                    z = "      " + str(z_cartesian_coords_list[conformer_number][i])
                else:
                    z = "     " + str(z_cartesian_coords_list[conformer_number][i])
                x += "0" * abs(18 - len(x))
                y += "0" * abs(18 - len(y))
                z += "0" * abs(18 - len(z))
                xyz_block += "\n  " + str(element_list[conformer_number][i]) + x + y + z
        mol = MolFromXYZBlock(xyz_block)
        return mol

    # Define function for finding optimal atom mapping for aligned conformers, using the Hungarian algorithm


    def maximum_atom_deviation_calculator(conf_a, conf_b):
        """Function to find the maximum atom deviation between two conformers.
        1) Aligns conformers by principal coordinates, based on their moments of inertia.
        2) Finds optimal atom mapping (lowest atom deviations) for two aligned conformers
        (considering all 8 axis reflections), using the Hungarian algorithm.
        3) Returns the maximum atom deviation for this optimal atom mapping."""
        # Align conformers by principal coordinates, based on their moments of inertia

        mol_a = RDKitMoleculeMaker(conf_number_a)
        mol_b = RDKitMoleculeMaker(conf_number_b)
        rdMolTransforms.CanonicalizeMol(mol_a)
        rdMolTransforms.CanonicalizeMol(mol_b)

        # Construct a list of lists of atom coordinates for conformer A - one list per element.

        unique_elements = set(atom.GetSymbol() for atom in mol_a.GetAtoms())
        mol_a_all_element_coords_list = []
        for element in unique_elements:
            element_coords_list = []
            for l, atom in enumerate(mol_a.GetAtoms()):
                position = mol_a.GetConformer().GetAtomPosition(l)
                if atom.GetSymbol() == element:
                    element_coords_list.append([position.x, position.y, position.z])
            mol_a_all_element_coords_list.append(element_coords_list)

        # Calculate atomic distances between aligned geometries using Hungarian algorithm, considering axis reflections

        all_mads = []
        for i in range(1, -2, -2):
            for j in range(1, -2, -2):
                for k in range(1, -2, -2):
                    aligned_conformer_atom_deviations = []
                    for x, element in enumerate(unique_elements):
                        mol_b_coords_list = []
                        for l, atom in enumerate(mol_b.GetAtoms()):
                            position = mol_b.GetConformer().GetAtomPosition(l)
                            if atom.GetSymbol() == element:
                                mol_b_coords_list.append([(position.x * i), (position.y * j), (position.z * k)])
                        mol_a_coords_array = np.array(mol_a_all_element_coords_list[x])
                        mol_b_coords_array = np.array(mol_b_coords_list)
                        distance_matrix = scipy_distance.cdist(mol_a_coords_array, mol_b_coords_array, 'euclidean')
                        row_ind, col_ind = linear_sum_assignment(distance_matrix)
                        atom_deviations = distance_matrix[row_ind, col_ind]
                        aligned_conformer_atom_deviations.extend(atom_deviations.tolist())
                    all_mads.append(max(aligned_conformer_atom_deviations))
        return min(all_mads)

    # Find pairs of conformers with similar or identical energies

    energy_threshold = energy_threshold / 627.5094740631  # Convert energy threshold from kcal/mol to hartrees
    number_of_confs = len(energies)
    if settings["Skip excluding duplicate conformers from input files made from output files"] is False:
        for conf_number_a, energy_a in enumerate(energies):
            removals = 0
            for conf_number_b in range(number_of_confs):
                conf_number_b = conf_number_b - removals
                if conf_number_b + 1 > len(energies):
                    break
                if conf_number_a < conf_number_b:
                    if abs(energy_a - energies[conf_number_b]) < energy_threshold:
                        maximum_atom_deviation = maximum_atom_deviation_calculator(conf_number_a, conf_number_b)
                        if maximum_atom_deviation < mad_threshold:
                            # Record which conformer(s) were duplicates

                            if list_of_conformer_suffixes:
                                dup_conf_name = list_of_conformer_suffixes[conf_number_b]
                                duplicate_conformers.append(dup_conf_name)
                                duplicated_conf_name = list_of_conformer_suffixes[conf_number_a]
                            else:
                                restart_loop = True  # Count the number of removed confs with lower conf numbers
                                new_dup_conf_number = conf_number_b + 1
                                y = []
                                while restart_loop is True:
                                    x = 0
                                    for number in dup_conf_numbers:
                                        if number not in y and number <= new_dup_conf_number:
                                            y.append(number)
                                            new_dup_conf_number += 1
                                            x += 1
                                            restart_loop = True
                                            break
                                    if x == 0:
                                        restart_loop = False
                                restart_loop = True  # Count the number of removed confs with lower conf numbers
                                new_duplicated_conf_number = conf_number_a + 1
                                y = []
                                while restart_loop is True:
                                    x = 0
                                    for number in dup_conf_numbers:
                                        if number not in y and number <= new_duplicated_conf_number:
                                            y.append(number)
                                            new_duplicated_conf_number += 1
                                            x += 1
                                            restart_loop = True
                                            break
                                    if x == 0:
                                        restart_loop = False
                                dup_conf_name = new_dup_conf_number
                                dup_conf_numbers.append(dup_conf_name)
                                duplicated_conf_name = new_duplicated_conf_number
                                if str(dup_conf_name).endswith("1") and not str(dup_conf_name).endswith("11"):
                                    dup_conf_name = str(dup_conf_name) + "st"
                                elif str(dup_conf_name).endswith("2") and not str(dup_conf_name).endswith("12"):
                                    dup_conf_name = str(dup_conf_name) + "nd"
                                elif str(dup_conf_name).endswith("3") and not str(dup_conf_name).endswith("13"):
                                    dup_conf_name = str(dup_conf_name) + "rd"
                                else:
                                    dup_conf_name = str(dup_conf_name) + "th"
                                if str(duplicated_conf_name).endswith("1") and not str(duplicated_conf_name).endswith(
                                    "11"
                                ):
                                    duplicated_conf_name = str(duplicated_conf_name) + "st"
                                elif str(duplicated_conf_name).endswith("2") and not str(duplicated_conf_name).endswith(
                                    "12"
                                ):
                                    duplicated_conf_name = str(duplicated_conf_name) + "nd"
                                elif str(duplicated_conf_name).endswith("3") and not str(duplicated_conf_name).endswith(
                                    "13"
                                ):
                                    duplicated_conf_name = str(duplicated_conf_name) + "rd"
                                else:
                                    duplicated_conf_name = str(duplicated_conf_name) + "th"
                                duplicate_conformers.append(dup_conf_name)
                            if settings["Duplicate conformer details"] is True:
                                if list_of_conformer_suffixes:
                                    if settings["Boltz energy type"] == "Gibbs free energy":
                                        dup_conf_text = (
                                            str(dup_conf_name)
                                            + " is a duplicate of "
                                            + str(duplicated_conf_name)
                                            + ":\nΔG = "
                                            + str((energies[conf_number_a] - energies[conf_number_b]) * 627.5094740631)
                                            + " kcal/mol\nMAD = "
                                            + str(maximum_atom_deviation)
                                            + " angstroms\n"
                                        )
                                    else:
                                        dup_conf_text = (
                                                str(dup_conf_name)
                                                + " is a duplicate of "
                                                + str(duplicated_conf_name)
                                                + ":\nΔE = "
                                                + str(
                                            (energies[conf_number_a] - energies[conf_number_b]) * 627.5094740631)
                                                + " kcal/mol\nMAD = "
                                                + str(maximum_atom_deviation)
                                                + " angstroms\n"
                                        )
                                else:
                                    if settings["Boltz energy type"] == "Gibbs free energy":
                                        dup_conf_text = (
                                            str(dup_conf_name)
                                            + " conformer is a duplicate of "
                                            + str(duplicated_conf_name)
                                            + " conformer:\nΔG = "
                                            + str((energies[conf_number_a] - energies[conf_number_b]) * 627.5094740631)
                                            + " kcal/mol\nMAD = "
                                            + str(maximum_atom_deviation)
                                            + " angstroms\n"
                                        )
                                    else:
                                        dup_conf_text = (
                                                str(dup_conf_name)
                                                + " conformer is a duplicate of "
                                                + str(duplicated_conf_name)
                                                + " conformer:\nΔE = "
                                                + str(
                                            (energies[conf_number_a] - energies[conf_number_b]) * 627.5094740631)
                                                + " kcal/mol\nMAD = "
                                                + str(maximum_atom_deviation)
                                                + " angstroms\n"
                                        )
                                duplicate_conf_details.append(dup_conf_text)
                            # Remove duplicate conformer from data

                            del energies[conf_number_b]
                            del element_list[conf_number_b]
                            del frequencies_list[conf_number_b]
                            del ir_intensities_list[conf_number_b]
                            del x_cartesian_coords_list[conf_number_b]
                            del y_cartesian_coords_list[conf_number_b]
                            del z_cartesian_coords_list[conf_number_b]
                            if list_of_conformer_suffixes:
                                del list_of_conformer_suffixes[conf_number_b]
                            if chk_conf_suffixes:
                                del chk_conf_suffixes[conf_number_b]
                            if wavelength_list:
                                del wavelength_list[conf_number_b]
                                del rotatory_strength_list[conf_number_b]
                                del oscillator_strength_list[conf_number_b]
                            if shielding_tensors:
                                del shielding_tensors[conf_number_b]
                            if frequency_rotatory_strengths:
                                del frequency_rotatory_strengths[conf_number_b]
                                del frequency_dipole_strengths[conf_number_b]
                            if optrot_wavelengths:
                                del optrot_strengths[conf_number_b]
                                del optrot_wavelengths[conf_number_b]
                            removals += 1  # This is because indices of remaining conformers are now reduced by 1,
                            # by removing the last conformer
    a = len(duplicate_conformers)
    if a > 0:
        b = ""
        c = ""
        for i, dup in enumerate(duplicate_conformers):
            dup = dup.removesuffix(".out")
            dup = dup.removesuffix(".log")
            b += dup
            if i + 1 != a:
                b += ", "
        if a > 1:
            c = "s"
        duplicate_conformers = [str(a), b, c]
    # Find relative energies of conformers in kJ/mol

    contributions_total = 0
    for conformer_number, conformer in enumerate(energies):
        if settings["Relative energy unit"] == "kcal/mol":
            relative_energy = (energies[conformer_number] - min(energies)) * hartree_to_kcalmol
        elif settings["Relative energy unit"] == "kJ/mol":
            relative_energy = (energies[conformer_number] - min(energies)) * hartree_to_kJmol
        if relative_energy == 0:
            relative_energy = int(relative_energy)
        relative_energies.insert(conformer_number, relative_energy)
        # Find Boltzmann-weighted population contributions of each conformer at room temperature

        population_contribution = exp(-1 * relative_energy / (boltz_constant * temperature))
        population_contributions.append(population_contribution)
        contributions_total += population_contribution
    # Calculate these population contributions as a percentage

    for contribution in population_contributions:
        weighting = contribution / contributions_total
        boltz_weights.append(weighting)
    # Calculate Boltzmann-averaged shielding tensors, if shielding tensors are present

    if shielding_tensors:
        # Calculate the shielding tensor component for each atom of each conformer

        for conformer in range(len(shielding_tensors)):
            components = []
            for shielding_tensor in shielding_tensors[conformer]:
                component = float(shielding_tensor) * boltz_weights[conformer]
                components.append(component)
            shielding_tensor_components.append(components)
        # Add together the shielding tensor components for each atom to produce Boltzmann-averaged shielding tensors

        for atom in range(len(shielding_tensor_components[0])):
            boltz_shielding_tensor = 0
            for conformer in range(len(shielding_tensor_components)):
                boltz_shielding_tensor += shielding_tensor_components[conformer][atom]
            boltz_shielding_tensors.append(boltz_shielding_tensor)
        # Pull out carbon and hydrogen shielding tensors into separate lists and calculate chemical shifts if selected

        row = []
        number_Cs = 0
        number_Hs = 0
        for atom in range(len(boltz_shielding_tensors)):
            if element_list[0][atom] == "C":
                number_Cs += 1
                label1 = label2 = "C"
                label1 += str(atom + 1)
                row.append(label1)
                label2 += str(number_Cs)
                row.append(label2)
                row.append(str(boltz_shielding_tensors[atom]))
                if c_slope and c_intercept is not None:
                    c_slope = float(c_slope)
                    c_intercept = float(c_intercept)
                    chemical_shift = (c_intercept - boltz_shielding_tensors[atom]) / -c_slope
                    row.append(str(chemical_shift))
                boltz_c_tensors.append(row)
            elif element_list[0][atom] == "H":
                number_Hs += 1
                label1 = label2 = "H"
                label1 += str(atom + 1)
                row.append(label1)
                label2 += str(number_Hs)
                row.append(label2)
                row.append(str(boltz_shielding_tensors[atom]))
                if h_slope and h_intercept is not None:
                    h_slope = float(h_slope)
                    h_intercept = float(h_intercept)
                    chemical_shift = (h_intercept - boltz_shielding_tensors[atom]) / -h_slope
                    row.append(str(chemical_shift))
                boltz_h_tensors.append(row)
            row = []
    # Calculate Boltzmann-averaged frequencies

    for conformer in range(len(frequencies_list)):
        components = []
        for frequency in frequencies_list[conformer]:
            component = float(frequency) * boltz_weights[conformer]
            components.append(component)
        frequency_components.append(components)
    # Add together the components for each frequency to produce Boltzmann-averaged frequencies

    for frequency in range(len(frequency_components[0])):
        boltz_frequency = 0
        for conformer in range(len(frequency_components)):
            boltz_frequency += frequency_components[conformer][frequency]
        boltz_frequencies.append(boltz_frequency)
    # Calculate Boltzmann-averaged IR intensities

    for conformer in range(len(ir_intensities_list)):
        components = []
        for ir_intensity in ir_intensities_list[conformer]:
            component = float(ir_intensity) * boltz_weights[conformer]
            components.append(component)
        ir_intensity_components.append(components)
    # Add together the components for each IR intensity to produce Boltzmann-averaged IR intensities

    for ir_intensity in range(len(ir_intensity_components[0])):
        boltz_ir_intensity = 0
        for conformer in range(len(ir_intensity_components)):
            boltz_ir_intensity += ir_intensity_components[conformer][ir_intensity]
        boltz_ir_intensities.append(boltz_ir_intensity)
    # If TD-DFT data present, calculate Boltzmann-averages

    if wavelength_list:
        # Calculate the wavelength component for each excited state of each conformer

        for conformer in range(len(wavelength_list)):
            components = []
            for excited_state in wavelength_list[conformer]:
                component = float(excited_state) * boltz_weights[conformer]
                components.append(component)
            wavelength_components.append(components)
        # Calculate the rotatory strength component for each excited state of each conformer

        for conformer in range(len(rotatory_strength_list)):
            components = []
            for excited_state in rotatory_strength_list[conformer]:
                component = float(excited_state) * boltz_weights[conformer]
                components.append(component)
            rotatory_strength_components.append(components)
        # Calculate the oscillator strength component for each excited state of each conformer

        for conformer in range(len(oscillator_strength_list)):
            components = []
            for excited_state in oscillator_strength_list[conformer]:
                component = float(excited_state) * boltz_weights[conformer]
                components.append(component)
            oscillator_strength_components.append(components)
        # Add together the components for each excited state to produce Boltzmann-averaged wavelengths

        for excited_state in range(len(wavelength_components[0])):
            boltz_wavelength = 0
            boltz_rotatory_strength = 0
            boltz_oscillator_strength = 0
            for conformer in range(len(wavelength_components)):
                boltz_wavelength += wavelength_components[conformer][excited_state]
                boltz_rotatory_strength += rotatory_strength_components[conformer][excited_state]
                boltz_oscillator_strength += oscillator_strength_components[conformer][excited_state]
            boltz_wavelengths.append(boltz_wavelength)
            boltz_rotatory_strengths.append(boltz_rotatory_strength)
            boltz_oscillator_strengths.append(boltz_oscillator_strength)
    # If VCD data present, calculate Boltzmann-averages

    if frequency_rotatory_strengths:
        # Calculate the (VCD) frequency rotatory strength component for each frequency of each conformer #

        for conformer in range(len(frequency_rotatory_strengths)):
            components = []
            for frequency in frequency_rotatory_strengths[conformer]:
                component = float(frequency) * boltz_weights[conformer]
                components.append(component)
            frequency_rotatory_strength_components.append(components)
        # Calculate the IR frequency dipole strength component for each frequency of each conformer

        for conformer in range(len(frequency_dipole_strengths)):
            components = []
            for frequency in frequency_dipole_strengths[conformer]:
                component = float(frequency) * boltz_weights[conformer]
                components.append(component)
            frequency_dipole_strength_components.append(components)
        # Add together the components for each frequency to produce Boltzmann-averaged frequencies

        for frequency in range(len(frequency_rotatory_strength_components[0])):
            boltz_frequency_rotatory_strength = 0
            boltz_frequency_dipole_strength = 0
            for conformer in range(len(frequency_components)):
                boltz_frequency_rotatory_strength += frequency_rotatory_strength_components[conformer][frequency]
                boltz_frequency_dipole_strength += frequency_dipole_strength_components[conformer][frequency]
            boltz_frequency_rotatory_strengths.append(boltz_frequency_rotatory_strength)
            boltz_frequency_dipole_strengths.append(boltz_frequency_dipole_strength)
    # If optical rotation data present, calculate Boltzmann-averages

    if optrot_wavelengths:
        # Calculate the optical rotation strength component for each excited state of each conformer
        #  consider wavelengths

        for conformer in range(len(optrot_strengths)):
            components = []
            for strength in optrot_strengths[conformer]:
                component = float(strength) * boltz_weights[conformer]
                components.append(component)
            optrot_strength_components.append(components)
        # Add together the components for each optical rotation strength to produce Boltzmann-averaged optical
        # rotation strengths

        for optrot_strength in range(len(optrot_strength_components[0])):
            boltz_optrot_strength = 0
            for conformer in range(len(optrot_strength_components)):
                boltz_optrot_strength += optrot_strength_components[conformer][optrot_strength]
            boltz_optrot_strengths.append(boltz_optrot_strength)
    # Find order of conformers ordered by ascending energy, dealing with spatially different conformers of
    # identical energies, if present

    energies_in_order = energies.copy()
    energies_in_order.sort()
    conformer_order = []
    previous_energy = 42  # I chose an arbitrary positive number here.
    for energy in energies_in_order:
        if energy != previous_energy:
            conformer_order.append(energies.index(energy))
            a = 1
        elif energy == previous_energy:
            start_at = energies.index(previous_energy) + a
            conformer_order.append(energies.index(energy, start_at))
            a += 1
        previous_energy = energy
    # Order the conformers in the data lists

    for conformer in conformer_order:
        ordered_energies.append(energies[conformer])
        ordered_relative_energies.append(relative_energies[conformer])
        ordered_boltz_weights.append(boltz_weights[conformer])
        ordered_population_contributions.append(population_contributions[conformer])
        ordered_element_list.append(element_list[conformer])
        ordered_x_cartesian_coords_list.append(x_cartesian_coords_list[conformer])
        ordered_y_cartesian_coords_list.append(y_cartesian_coords_list[conformer])
        ordered_z_cartesian_coords_list.append(z_cartesian_coords_list[conformer])
        ordered_frequencies_list.append(frequencies_list[conformer])
        ordered_ir_intensities_list.append(ir_intensities_list[conformer])
        if wavelength_list:
            ordered_wavelength_list.append(wavelength_list[conformer])
            ordered_rotatory_strengths.append(rotatory_strength_list[conformer])
            ordered_oscillator_strengths.append(oscillator_strength_list[conformer])
        if shielding_tensors:
            ordered_shielding_tensors.append(shielding_tensors[conformer])
        if frequency_rotatory_strengths:
            ordered_frequency_rotatory_strengths_list.append(frequency_rotatory_strengths[conformer])
            ordered_frequency_dipole_strengths_list.append(frequency_dipole_strengths[conformer])
        if optrot_wavelengths:
            ordered_optrot_list.append(optrot_strengths[conformer])
            ordered_optrot_wavelengths_list.append(optrot_wavelengths[conformer])
    return (
        ordered_energies,
        ordered_relative_energies,
        ordered_boltz_weights,
        ordered_element_list,
        ordered_x_cartesian_coords_list,
        ordered_y_cartesian_coords_list,
        ordered_z_cartesian_coords_list,
        ordered_wavelength_list,
        ordered_rotatory_strengths,
        ordered_shielding_tensors,
        boltz_shielding_tensors,
        boltz_c_tensors,
        boltz_h_tensors,
        boltz_wavelengths,
        boltz_rotatory_strengths,
        sp_functional_and_basis_set,
        sp_solvent,
        sp_dispersion,
        opt_freq_functional_and_basis_set,
        opt_freq_solvent,
        opt_freq_dispersion,
        nmr_functional_and_basis_set,
        nmr_solvent,
        nmr_dispersion,
        tddft_functional_and_basis_set,
        tddft_solvent,
        tddft_dispersion,
        results_directory,
        boltz_oscillator_strengths,
        boltz_frequencies,
        ordered_frequencies_list,
        ordered_oscillator_strengths,
        duplicate_conformers,
        data_analysis_error_check,
        ordered_frequency_rotatory_strengths_list,
        boltz_frequency_rotatory_strengths,
        ordered_optrot_wavelengths_list,
        ordered_optrot_list,
        boltz_optrot_strengths,
        ordered_frequency_dipole_strengths_list,
        boltz_frequency_dipole_strengths,
        or_functional_and_basis_set,
        or_solvent,
        or_dispersion,
        duplicate_conf_details,
        relative_energies,
        element_list,
        x_cartesian_coords_list,
        y_cartesian_coords_list,
        z_cartesian_coords_list,
        list_of_conformer_suffixes,
        boltz_ir_intensities,
        ordered_ir_intensities_list,
        number_imag_freq_confs_removed,
        chk_conf_suffixes,
    )
