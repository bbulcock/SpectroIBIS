# Module containing functions for parsing, extracting and checking a) computational data from Gaussian/ORCA output files,
# and b) conformer Cartesian coordinates from .xyz/.sdf files.


from re import search, findall, sub, IGNORECASE, DOTALL, MULTILINE
from rdkit.Chem import MolFromXYZBlock, rdMolTransforms  # Version 2025.9.3
import numpy as np  # Version 2.2.6
from scipy.optimize import linear_sum_assignment  # Version 1.15.3
from scipy.spatial import distance as scipy_distance


# Parses Gaussian or ORCA output files. Uses a list of filenames as input.


def parse(list_of_filepaths, settings):
    """Extracts key data from a Gaussian/ORCA output files (as a list of file paths),
    then performs a variety of error checks."""
    # Define lists

    energies = []
    element_list = []
    x_cartesian_coords_list = []
    y_cartesian_coords_list = []
    z_cartesian_coords_list = []
    sp_energies = []
    opt_freq_energies = []
    wavelength_list = []
    rotatory_strength_list = []
    oscillator_strength_list = []
    shielding_tensors = []
    frequencies = []
    ir_intensities = []
    frequency_rotatory_strengths = []
    frequency_dipole_strengths = []
    optrot_strengths = []
    optrot_wavelengths = []
    Gibbs_corrections = []
    conformer_element_list = []
    conformer_x_cartesian_coords_list = []
    conformer_y_cartesian_coords_list = []
    conformer_z_cartesian_coords_list = []
    all_coordinates_texts = []
    list_of_conformer_suffixes = []
    min_frequencies = []
    list_of_file_contents = []
    imaginary_freq_confs = []
    chk_conf_suffixes = []
    all_sp_calc_details = []
    all_opt_calc_details = []
    all_nmr_calc_details = []
    all_or_calc_details = []
    all_tddft_calc_details = []

    # Define variables

    sp_functional_and_basis_set = ""
    sp_solvent = ""
    sp_dispersion = ""
    opt_functional_and_basis_set = ""
    opt_solvent = ""
    opt_dispersion = ""
    nmr_functional_and_basis_set = ""
    nmr_solvent = ""
    nmr_dispersion = ""
    tddft_functional_and_basis_set = ""
    tddft_solvent = ""
    tddft_dispersion = ""
    or_functional_and_basis_set = ""
    or_solvent = ""
    or_dispersion = ""
    e = ""
    imaginary_freq_confs_text = ""
    are_there_files_without_conf_suffix = ""
    calc_software = ""
    # Define regex patterns

    gaussian_calc_details_functional_and_basis_set_regex = r"\S+/\S+"
    gaussian_calc_details_solvent_regex = r"solvent=(\w+-*\w*)"
    gaussian_calc_details_dispersion_regex = r"EmpiricalDispersion=(\w+)"
    gaussian_sp_calc_details_regex = r"---\n( #\w? .{0,200} SP (.){0,200})\n ---"
    gaussian_opt_calc_details_regex = r"---\n( #\w? .{0,200}opt(.){0,200})\n ---"
    gaussian_nmr_calc_details_regex = r"---\n( #\w? .{0,200}NMR(.){0,200})\n ---"
    gaussian_or_calc_details_regex = r"---\n( #\w? .{0,200}Polar=OptRot(.){0,200})\n ---"
    gaussian_tddft_calc_details_regex = r"---\n( #\w? .{0,200}td=?\(nstates?=(.){0,200})\n ---"
    gaussian_coords_block_regex = r"1\\1\\\S*\\FOpt\\.*?\\\\@|1\|1\|\S*\|FOpt\|.*?\|\|@"
    gaussian_coords_regex = r"([A-Z][a-z]?),(-?\d+.\d+),(-?\d+.\d+),(-?\d+.\d+)"
    gaussian_gibbs_free_energies_regex = r"Sum of electronic and thermal Free Energies=\s*(-[0-9]+\.[0-9]+).*"
    gaussian_sp_energies_regex = r"=(-\d+\.\d+)\\RMSD=\d\.\d+e"
    gaussian_freq_section_regex = r"Frequencies --.*?- Thermochemistry -"
    gaussian_frequencies_regex = r"Frequencies --\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
    gaussian_ir_intensities_regex = r"IR Inten    --\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
    gaussian_frequency_rotatory_strengths_regex = r"Rot\. str\.   --\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
    gaussian_frequency_dipole_strengths_regex = r"Dip\. str\.   --\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
    gaussian_nmr_regex = r" Calculating GIAO nuclear magnetic shielding tensors\..[^*]*\*\*\*\*\*\*\*\*\*\*"
    gaussian_shielding_tensors_regex = r".*Isotropic =\s*(-?[0-9]+\.[0-9]+).*"
    gaussian_ecd_regex = (
        r" R\(length\).*? SavETr:  write IOETrn=   \d+ NScale= \d+ NData=  \d+ NLR=1 "
    )
    gaussian_rotatory_strength_regex = r" +\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +(-?\d+\.\d+)\n"
    gaussian_optrot_section_regex = r"\[Alpha\] \( \d+\.\d+ A\) = +-?\d+\.\d+ deg\..[^*]*\*\*\*\*\*\*\*\*\*\*"
    gaussian_oscillator_strength_regex = (
        r"Excited State +\d+:      Singlet-A      \d+\.\d+ eV  \d+\.\d\d nm  f=(\d+\.\d+)  <S\*\*2>=0\.000"
    )
    gaussian_wavelength_regex = (
        r"Excited State +\d+:      Singlet-A      \d+\.\d+ eV  (\d+\.\d\d) nm  f=\d+\.\d+  <S\*\*2>=0\.000"
    )
    gaussian_optrot_regex = r"\[Alpha\] \( (\d+\.\d+) A\) = +(-?\d+\.\d+) deg\."
    gaussian_Gibbs_corrections_regex = r"\n Thermal correction to Gibbs Free Energy= *(\d+\.\d+)"

    orca_calc_details_functional_and_basis_set_regex = r"^ *(\S+ +\S+)"
    orca_calc_details_solvent_regex = r"CPCMC?\((\w+-*\w*)\)"
    orca_calc_details_dispersion_regex = r"D4|D3BJ|D3ZERO|D2|NL|SCNL"
    orca_calc_details_regex = r"(\n\| *\d+> +!.[^*]*>) +\*"
    orca_sp_calc_details_regex = r"(\n\| *\d+> +!.[^*]*( SP | energy ).[^*]*>) +(\*)"
    orca_opt_calc_details_regex = r"(\n\| *\d+> +!.[^*]*Opt.[^*]*>) +(\*)"
    orca_nmr_calc_details_regex = r"(\n\| *\d+> +!.[^*]*NMR.[^*]*>) +(\*)"
    orca_or_calc_details_regex = r""  # optical rotation not available in ORCA 6
    orca_tddft_calc_details_regex = r"(\n\| *\d+> +!.[^*]*%TDDFT.[^*]*>) +(\*)"
    orca_sp_energies_regex = r"----\nFINAL SINGLE POINT ENERGY +(-\d+\.\d+)"
    orca_freq_section_regex = r"VIBRATIONAL FREQUENCIES\n------.*?NORMAL MODES\n---"
    orca_frequencies_regex = r" *\d+: *(-?\d+\.\d+) cm\*\*-1"
    orca_ir_section_regex = r"IR SPECTRUM\n------.*?---\nTHERMOCHEMISTRY AT"
    orca_ir_blocks_regex = (
        r"\n *\d+: +\d+\.\d+ +\d+\.\d+ +(\d+\.\d+) +\d+\.\d+ +\(-? ?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+\)"
    )
    orca_ir_intensities_regex = r" *\d+: *(-?\d+\.\d+) km/mol"
    orca_vcd_section_regex = r"VCD SPECTRUM CALCULATION\n------.*?Maximum memory used throughout the entire "
    orca_frequency_rotatory_strengths_regex = r"\n *\d+ +\d+\.\d+ +(-?\d+\.\d+)"
    orca_frequency_extinction_coefficients_regex = (
        r"\n *\d+: +\d+\.\d+ +(\d+\.\d+) +\d+\.\d+ +\d+\.\d+ +\(-? ?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+\)"
    )
    orca_coords_block_regex = (
        r"\*\*\* FINAL ENERGY EVALUATION AT THE STATIONARY POINT \*\*\*.*?CARTESIAN COORDINATES \(A\.U\.\)"
    )
    orca_coords_regex = r"([A-Z][a-z]?)\s+(-?\d+.\d+)\s+(-?\d+.\d+)\s+(-?\d+.\d+)"
    orca_gibbs_free_energies_regex = r"\nFinal Gibbs free energy *\.\.\. *(-\d+\.\d+) Eh"
    orca_nmr_regex = (
        r"CHEMICAL SHIELDING SUMMARY \(ppm\).*?NMR"
    )
    orca_shielding_tensors_regex = r"\n\s+\d+\s+\w+\s+(-?\d+\.\d+)+"
    orca_ecd_regex = r"CD SPECTRUM\n-------------------.*?\n\n"
    orca_rotatory_strength_regex = r"\n\s+\d+\s+-?\d+\.\d+\s+\d+\.\d+\s+(-?\d+\.\d+)"
    orca_oscillator_strength_regex = (
        r"\n\s+\d+\s+\d+\.\d+\s+\d+\.\d+\s+(-?\d+\.\d+)\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+"
    )
    orca_wavelength_regex = r"\n\s+\d+\s+-?\d+\.\d+\s+(\d+\.\d+)"
    orca_optrot_regex = r""  # optical rotation not available in ORCA 6.1
    orca_optrot_section_regex = r""  # optical rotation not available in ORCA 6.1
    orca_Gibbs_corrections_regex = r"\nG-E\(el\) +\.\.\. *(\d+\.\d+) Eh"

    # Set default value for error status

    parser_error_check = "No error detected"
    error_message = ""

    # Check selected files are all *.out or *.log files

    for filename in list_of_filepaths:
        if filename.endswith(".out") is False and filename.endswith(".log") is False:
            parser_error_check = "Error detected"
            error_message = "Incorrect file format, or mixed file formats (e.g. *.out and *.xyz)."
            return (
                "",
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                parser_error_check,
                0,
                error_message,
            )
    # Produce a list of filenames ordered by conformer, if separate comp chem output files for the same conformer are
    # present.

    list_of_filepaths_by_conformer = []
    blacklist = []
    stop = False
    for index, filename in enumerate(list_of_filepaths):
        if stop == False:
            # Check for filename in blacklist

            ignore = False
            if blacklist:
                for blacklist_filename in blacklist:
                    if filename == blacklist_filename:
                        ignore = True
            if ignore == False:
                # Get filenames with same conformer suffix and group them by conformer number in their filename.

                conformer_suffix = search(
                    "conf-\d+\.log|conformer-\d+\.log|M\d\d\d\d\.log|conf-\d+\.out|conformer-\d+\.out|M\d\d\d\d\.out",
                    filename,
                    IGNORECASE
                )
                if conformer_suffix is not None:
                    conformer_suffix = conformer_suffix[0]
                    matches = []
                    for filename in list_of_filepaths:
                        match = findall(".*" + conformer_suffix, filename, DOTALL)
                        if match:
                            match = match[0]
                            matches.append(match)
                            blacklist.append(match)
                    list_of_filepaths_by_conformer.append(matches)
                if conformer_suffix is None:
                    are_there_files_without_conf_suffix = True
                    if len(list_of_filepaths) > 1:
                        list_of_filepaths_by_conformer.append(list_of_filepaths[index])
                    if len(list_of_filepaths) == 1:
                        list_of_filepaths_by_conformer.append(filename)
    # Open files and extract text.

    for filename_list_index, conformer in enumerate(list_of_filepaths_by_conformer):
        data = ""
        if isinstance(conformer, list):
            for file in conformer:
                with open(file) as file:
                    text = file.read()
                    data += text
            # Record conformer suffix numbers

            conformer_suffix = search(
                "(conf-\d+)\.log|(conformer-\d+)\.log|(M\d\d\d\d)\.log|(conf-\d+)\.out|(conformer-\d+)\.out|("
                "M\d\d\d\d)\.out",
                conformer[0],
                IGNORECASE
            )
            if conformer_suffix is not None:
                list_of_conformer_suffixes.append(conformer_suffix[0])
            elif conformer_suffix is None and len(list_of_filepaths_by_conformer) == 1:
                for directory in conformer:
                    name = directory.split("/")[-1]
                    list_of_conformer_suffixes.append(name)
        elif isinstance(conformer, str):
            with open(conformer) as file:
                data = file.read()
        # Distinguish if files are Gaussian or 0RCA, then assign relevant regex patterns.

        gaussian = False
        orca = False
        if search(r"\* O   R   C   A \*", data) is not None:
            orca = True
            calc_software = "orca"
            calc_details_functional_and_basis_set_regex = orca_calc_details_functional_and_basis_set_regex
            calc_details_solvent_regex = orca_calc_details_solvent_regex
            calc_details_dispersion_regex = orca_calc_details_dispersion_regex
            sp_calc_details_regex = orca_sp_calc_details_regex
            opt_calc_details_regex = orca_opt_calc_details_regex
            nmr_calc_details_regex = orca_nmr_calc_details_regex
            or_calc_details_regex = orca_or_calc_details_regex
            tddft_calc_details_regex = orca_tddft_calc_details_regex
            coords_block_regex = orca_coords_block_regex
            coords_regex = orca_coords_regex
            gibbs_free_energies_regex = orca_gibbs_free_energies_regex
            sp_energies_regex = orca_sp_energies_regex
            freq_section_regex = orca_freq_section_regex
            frequencies_regex = orca_frequencies_regex
            ir_intensities_regex = orca_ir_intensities_regex
            frequency_rotatory_strengths_regex = orca_frequency_rotatory_strengths_regex
            frequency_dipole_strengths_regex = orca_frequency_extinction_coefficients_regex
            # Note that while Gaussian computes IR dipole strengths, ORCA 6 actually computes IR extinction
            # coefficients, NOT dipole strengths. This difference is accounted for by SpectroIBIS later in the final
            # .ir.bil and .docx files.

            nmr_regex = orca_nmr_regex
            shielding_tensors_regex = orca_shielding_tensors_regex
            ecd_regex = orca_ecd_regex
            rotatory_strength_regex = orca_rotatory_strength_regex
            oscillator_strength_regex = orca_oscillator_strength_regex
            wavelength_regex = orca_wavelength_regex
            optrot_section_regex = orca_optrot_section_regex
            optrot_regex = orca_optrot_regex
            Gibbs_corrections_regex = orca_Gibbs_corrections_regex
        elif search(r"Gaussian, Inc\.  All Rights Reserved\.", data) is not None:
            gaussian = True
            calc_software = "gaussian"
            calc_details_functional_and_basis_set_regex = gaussian_calc_details_functional_and_basis_set_regex
            calc_details_solvent_regex = gaussian_calc_details_solvent_regex
            calc_details_dispersion_regex = gaussian_calc_details_dispersion_regex
            sp_calc_details_regex = gaussian_sp_calc_details_regex
            opt_calc_details_regex = gaussian_opt_calc_details_regex
            nmr_calc_details_regex = gaussian_nmr_calc_details_regex
            or_calc_details_regex = gaussian_or_calc_details_regex
            tddft_calc_details_regex = gaussian_tddft_calc_details_regex
            coords_block_regex = gaussian_coords_block_regex
            coords_regex = gaussian_coords_regex
            gibbs_free_energies_regex = gaussian_gibbs_free_energies_regex
            sp_energies_regex = gaussian_sp_energies_regex
            freq_section_regex = gaussian_freq_section_regex
            frequencies_regex = gaussian_frequencies_regex
            ir_intensities_regex = gaussian_ir_intensities_regex
            frequency_rotatory_strengths_regex = gaussian_frequency_rotatory_strengths_regex
            frequency_dipole_strengths_regex = gaussian_frequency_dipole_strengths_regex
            nmr_regex = gaussian_nmr_regex
            shielding_tensors_regex = gaussian_shielding_tensors_regex
            ecd_regex = gaussian_ecd_regex
            rotatory_strength_regex = gaussian_rotatory_strength_regex
            oscillator_strength_regex = gaussian_oscillator_strength_regex
            wavelength_regex = gaussian_wavelength_regex
            optrot_section_regex = gaussian_optrot_section_regex
            optrot_regex = gaussian_optrot_regex
            Gibbs_corrections_regex = gaussian_Gibbs_corrections_regex
        else:
            # Some ORCA 6 output files use UTF-16 encoding - will now try opening with this
            if isinstance(conformer, list):
                for file in conformer:
                    with open(file, encoding="utf-16") as file:
                        text = file.read()
                        data += text
                # Record conformer suffix numbers

                conformer_suffix = search(
                    "(conf-\d+)\.log|(conformer-\d+)\.log|(M\d\d\d\d)\.log|(conf-\d+)\.out|(conformer-\d+)\.out|("
                    "M\d\d\d\d)\.out",
                    conformer[0],
                    IGNORECASE
                )
                if conformer_suffix is not None:
                    list_of_conformer_suffixes.append(conformer_suffix[0])
                elif conformer_suffix is None and len(list_of_filepaths_by_conformer) == 1:
                    for directory in conformer:
                        name = directory.split("/")[-1]
                        list_of_conformer_suffixes.append(name)
            elif isinstance(conformer, str):
                with open(conformer, encoding="utf-16") as file:
                    data = file.read()

            if search(r"\* O   R   C   A \*", data) is not None:
                orca = True
                calc_software = "orca"
                calc_details_functional_and_basis_set_regex = orca_calc_details_functional_and_basis_set_regex
                calc_details_solvent_regex = orca_calc_details_solvent_regex
                calc_details_dispersion_regex = orca_calc_details_dispersion_regex
                sp_calc_details_regex = orca_sp_calc_details_regex
                opt_calc_details_regex = orca_opt_calc_details_regex
                nmr_calc_details_regex = orca_nmr_calc_details_regex
                or_calc_details_regex = orca_or_calc_details_regex
                tddft_calc_details_regex = orca_tddft_calc_details_regex
                coords_block_regex = orca_coords_block_regex
                coords_regex = orca_coords_regex
                gibbs_free_energies_regex = orca_gibbs_free_energies_regex
                sp_energies_regex = orca_sp_energies_regex
                freq_section_regex = orca_freq_section_regex
                frequencies_regex = orca_frequencies_regex
                ir_intensities_regex = orca_ir_intensities_regex
                frequency_rotatory_strengths_regex = orca_frequency_rotatory_strengths_regex
                frequency_dipole_strengths_regex = orca_frequency_extinction_coefficients_regex
                nmr_regex = orca_nmr_regex
                shielding_tensors_regex = orca_shielding_tensors_regex
                ecd_regex = orca_ecd_regex
                rotatory_strength_regex = orca_rotatory_strength_regex
                oscillator_strength_regex = orca_oscillator_strength_regex
                wavelength_regex = orca_wavelength_regex
                optrot_section_regex = orca_optrot_section_regex
                optrot_regex = orca_optrot_regex
                Gibbs_corrections_regex = orca_Gibbs_corrections_regex
            else:  # Could not recognise output file
                parser_error_check = "Error detected"
                error_message = (
                    "Could not recognise contents of output file(s).\nCheck file contents are standard for "
                    "Gaussian/ORCA output files. "
                )
                return (
                    "",
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    parser_error_check,
                    0,
                    error_message,
                )
        # Check if dedicated single-point energy calcs are present, extract theory level and energies

        sp_calc_details = findall(sp_calc_details_regex, data, IGNORECASE | DOTALL)
        if sp_calc_details:

            # Remove filenames from parsed calc details. Then, if they no longer contain keyword, remove them.
            if orca:
                for i in sp_calc_details:
                    unwanted_name = findall(r"%base \".*\"", i[0])
                    if unwanted_name:
                        x = i[0].replace(unwanted_name[0], "")
                        if ' sp ' not in x.lower() and ' energy ' not in x.lower():
                            sp_calc_details.remove(i)
            for i in sp_calc_details:
                all_sp_calc_details.append(i)
            if gaussian:
                # Replace sp_calc_details tuple with new tuple with cleaned data

                c = []
                for i, j in enumerate(sp_calc_details):
                    b = sp_calc_details[i][0].replace("\n ", "")
                    a = (b, sp_calc_details[i][1])
                    c.append(a)
                sp_calc_details = c
            sp_solvent = findall(calc_details_solvent_regex, sp_calc_details[0][0], IGNORECASE)
            if not sp_solvent and orca:
                sp_solvent = findall(r'SMDSOLVENT "(\w+-*\w*)"', sp_calc_details[0][0], IGNORECASE)
            if not sp_solvent:
                sp_solvent = ""
            elif sp_solvent:
                sp_solvent = sp_solvent[0]
            if orca:
                keyword_sp_calc_details = []
                for number, section in enumerate(sp_calc_details):
                    b = sub(r"\n\| *\d+> +[^!]*", r"", sp_calc_details[number][0])
                    b = sub("!", "", b)
                    b = sub("NOAUTOSTART|PAL\d+|xyzfile", "", b, IGNORECASE)

                    # Clean up ORCA 6's keywords for Martin's DSD functionals
                    b = sub("REVDSD-PBEP86", "revDSD-PBEP86", b, IGNORECASE)
                    b = sub("REVDOD-PBEP86", "revDOD-PBEP86", b, IGNORECASE)
                    b = sub("DSD-PBEP86/2013", "DSD-PBEP86", b, IGNORECASE)
                    b = sub("DSD-BLYP/2013", "DSD-BLYP", b, IGNORECASE)
                    b = sub("D-PBEP86-D4/2021", "D-PBEP86 D4", b, IGNORECASE)
                    b = sub("D-PBEP86/2021", "D-PBEP86", b, IGNORECASE)
                    keyword_sp_calc_details.append(b)
                sp_calc_details = [keyword_sp_calc_details]
            sp_dispersion = findall(calc_details_dispersion_regex, sp_calc_details[0][0], IGNORECASE)
            if not sp_dispersion:
                sp_dispersion = ""
            elif sp_dispersion:
                sp_dispersion = sp_dispersion[0]
            if gaussian:
                sp_calc_details = sub(
                    "\n ", "", sp_calc_details[0][0]
                )  # Remove "\n " in case text of interest is split across multiple lines
            if orca:
                sp_calc_details = sub(sp_dispersion, "", sp_calc_details[0][0])
            sp_functional_and_basis_set = findall(calc_details_functional_and_basis_set_regex, sp_calc_details)[0]
            if orca:
                sp_functional_and_basis_set = sub(r" +", r"/", sp_functional_and_basis_set)
            sp_functional_and_basis_set = sub(
                "#\w? *|opt|EmpiricalDispersion=\w*$|scrf=\w*$|Geom=\w*$|Guess=\w*$",
                "",
                sp_functional_and_basis_set,
                IGNORECASE,
            )
            #  Reformat functional and basis set

            sp_functional_and_basis_set = sp_functional_and_basis_set.replace("b3lyp", "B3LYP")
            sp_functional_and_basis_set = sp_functional_and_basis_set.replace("wB97XD", "ωB97X-D")
            sp_functional_and_basis_set = sp_functional_and_basis_set.replace("wb97xd", "ωB97X-D")
            sp_functional_and_basis_set = sp_functional_and_basis_set.replace("wB97MV", "ωB97M-V")
            sp_functional_and_basis_set = sp_functional_and_basis_set.replace("wb97mv", "ωB97M-V")
            sp_functional_and_basis_set = sp_functional_and_basis_set.replace("wB", "ωB")
            sp_functional_and_basis_set = sp_functional_and_basis_set.replace("M062X", "M06-2X")
            sp_functional_and_basis_set = sp_functional_and_basis_set.replace("m062x", "M06-2X")
            sp_functional_and_basis_set = sp_functional_and_basis_set.replace("Def2", "def2")
            if "def2" in sp_functional_and_basis_set and "def2-" not in sp_functional_and_basis_set:
                sp_functional_and_basis_set = sp_functional_and_basis_set.replace("def2", "def2-")
            sp_functional_and_basis_set = sp_functional_and_basis_set.replace("AUG-cc-", "aug-cc-")
            if sp_dispersion.casefold() == "gd2":
                sp_dispersion = "D2"
            elif sp_dispersion.casefold() == "gd3":
                sp_dispersion = "D3"
            elif sp_dispersion.casefold() == "gd3bj":
                sp_dispersion = "D3(BJ)"
            if sp_dispersion:
                sp_functional = sp_functional_and_basis_set.split("/")[0]
                sp_basis_set = sp_functional_and_basis_set.split("/")[1]
                sp_functional += "-" + sp_dispersion
                sp_functional_and_basis_set = sp_functional + "/" + sp_basis_set
            sp_calc_details = sp_calc_details.lstrip()
            # Find blocks of text containing sp energies corresponding to the dedicated sp calcs

            data_no_EOL = data.replace("\n ", "")
            sp_blocks = []
            if gaussian:
                end_sections = findall(r"1\\1\\.*?\\\\@|1\|1\|.*?\|\|@", data_no_EOL, DOTALL)
                for section in end_sections:
                    section = section.replace("\n ", "")
                    section = section.replace("|", "\\")
                    if sp_calc_details in section:
                        sp_blocks.append(section)
            elif orca:
                # Find job blocks

                job_blocks = findall(
                    r"---------------------\nBASIS SET INFORMATION\n---------------------.*?\nTimings for individual modules:", data, DOTALL
                )
                # Find job numbers corresponding to dedicated sp energy calcs

                calc_blocks = findall(orca_calc_details_regex, data, DOTALL)
                sp_calc_indices = []
                for index, section in enumerate(calc_blocks):
                    if " sp " in section.casefold() or "energy" in section.casefold():
                        sp_calc_indices.append(index)
                # Extract job blocks corresponding to dedicated sp calcs

                for index in sp_calc_indices:
                    sp_blocks.append(job_blocks[index])
            # Extract sp energies

            for sp_block in sp_blocks:
                sp_energy = findall(sp_energies_regex, sp_block)
                if len(sp_energy) == 1:  # Check that this job section has only one sp energy
                    sp_energies.append(float(sp_energy[0]))
        if settings["Boltz energy type"] == "Gibbs free energy":
            # Extract energy corrections and use these to calculate Gibbs free energies

            Gibbs_corrections_in_file = findall(Gibbs_corrections_regex, data)
            if Gibbs_corrections_in_file:
                for correction in Gibbs_corrections_in_file:
                    Gibbs_corrections.append(correction)
        # Find opt freq calc DFT theory level and solvation.

        if orca:
            data = data.replace("> %coords", "> *")  # For orca files, convert '> %coords' to '> *' for later regex
            # analysis.
        opt_calc_details = findall(opt_calc_details_regex, data, IGNORECASE | DOTALL)
        if opt_calc_details:

            # Remove filenames from parsed calc details. Then, if they no longer contain keyword, remove them.
            if orca:
                for i in opt_calc_details:
                        unwanted_name = findall(r"%base \".*\"", i[0])
                        if unwanted_name:
                            x = i[0].replace(unwanted_name[0], "")
                            if 'opt' not in x.lower():
                                opt_calc_details.remove(i)
            # Remove Polar=OptRot from parsed opt calc details. Then, if they no longer contain keyword, remove them.
            for i in opt_calc_details:
                x = i[0].lower().replace("polar=optrot", "")
                if 'opt' not in x:
                    opt_calc_details.remove(i)
            for i in opt_calc_details:
                all_opt_calc_details.append(i)
            if gaussian:

                # Replace opt_calc_details tuple with new tuple with cleaned data

                c = []
                for i, j in enumerate(opt_calc_details):
                    b = opt_calc_details[i][0].replace("\n ", "")
                    a = (b, opt_calc_details[i][1])
                    c.append(a)
                opt_calc_details = c
            opt_solvent = findall(calc_details_solvent_regex, opt_calc_details[0][0], IGNORECASE)
            if not opt_solvent and orca:
                opt_solvent = findall(r'SMDSOLVENT "(\w+-*\w*)"', opt_calc_details[0][0], IGNORECASE)
            if not opt_solvent:
                opt_solvent = ""
            elif opt_solvent:
                opt_solvent = opt_solvent[0]
            if orca:
                keyword_opt_calc_details = []
                for number, section in enumerate(opt_calc_details):
                    b = sub(r"\n\| *\d+> +[^!]*", r"", opt_calc_details[number][0])
                    b = sub("!", "", b)
                    b = sub("NOAUTOSTART|PAL\d+|xyzfile", "", b, IGNORECASE)

                    # Clean up ORCA 6's keywords for Martin's DSD functionals
                    b = sub("REVDSD-PBEP86", "revDSD-PBEP86", b, IGNORECASE)
                    b = sub("REVDOD-PBEP86", "revDOD-PBEP86", b, IGNORECASE)
                    b = sub("DSD-PBEP86/2013", "DSD-PBEP86", b, IGNORECASE)
                    b = sub("DSD-BLYP/2013", "DSD-BLYP", b, IGNORECASE)
                    b = sub("D-PBEP86-D4/2021", "D-PBEP86 D4", b, IGNORECASE)
                    b = sub("D-PBEP86/2021", "D-PBEP86", b, IGNORECASE)
                    keyword_opt_calc_details.append(b)
                opt_calc_details = [keyword_opt_calc_details]
            opt_dispersion = findall(calc_details_dispersion_regex, opt_calc_details[0][0], IGNORECASE)
            if not opt_dispersion:
                opt_dispersion = ""
            elif opt_dispersion:
                opt_dispersion = opt_dispersion[0]
            if gaussian:
                opt_calc_details = sub("\n ", "", opt_calc_details[0][0])  # Remove "\n " in case text of interest is
                # split across multiple lines
            if orca:
                opt_calc_details = sub(opt_dispersion, "", opt_calc_details[0][0])
            opt_functional_and_basis_set = findall(calc_details_functional_and_basis_set_regex, opt_calc_details)[0]
            if orca:
                opt_functional_and_basis_set = sub(r" +", r"/", opt_functional_and_basis_set)
                opt_functional = opt_functional_and_basis_set.lstrip().split("/")
                opt_functional = opt_functional[0]
                if opt_functional.casefold().endswith(
                    "-3c"
                ):  # -3c methods don't have a basis set, so remove the next word that would otherwise be assumed to be a basis set.
                    opt_functional_and_basis_set = opt_functional
            opt_functional_and_basis_set = sub(
                "#\w? *|opt|EmpiricalDispersion=\w*$|scrf=\w*$|Geom=\w*$|Guess=\w*$",
                "",
                opt_functional_and_basis_set,
                IGNORECASE,
            )
            #  Reformat functional and basis set

            opt_functional_and_basis_set = opt_functional_and_basis_set.replace("b3lyp", "B3LYP")
            opt_functional_and_basis_set = opt_functional_and_basis_set.replace("wB97XD", "ωB97X-D")
            opt_functional_and_basis_set = opt_functional_and_basis_set.replace("wb97xd", "ωB97X-D")
            opt_functional_and_basis_set = opt_functional_and_basis_set.replace("wB97MV", "ωB97M-V")
            opt_functional_and_basis_set = opt_functional_and_basis_set.replace("wb97mv", "ωB97M-V")
            opt_functional_and_basis_set = opt_functional_and_basis_set.replace("wB", "ωB")
            opt_functional_and_basis_set = opt_functional_and_basis_set.replace("M062X", "M06-2X")
            opt_functional_and_basis_set = opt_functional_and_basis_set.replace("m062x", "M06-2X")
            opt_functional_and_basis_set = opt_functional_and_basis_set.replace("Def2", "def2")
            if "def2" in opt_functional_and_basis_set and "def2-" not in opt_functional_and_basis_set:
                opt_functional_and_basis_set = opt_functional_and_basis_set.replace("def2", "def2-")
            opt_functional_and_basis_set = opt_functional_and_basis_set.replace("AUG-cc-", "aug-cc-")
            if opt_dispersion.casefold() == "gd2":
                opt_dispersion = "D2"
            elif opt_dispersion.casefold() == "gd3":
                opt_dispersion = "D3"
            elif opt_dispersion.casefold() == "gd3bj":
                opt_dispersion = "D3(BJ)"
            if opt_dispersion:
                opt_functional = opt_functional_and_basis_set.split("/")[0]
                opt_basis_set = opt_functional_and_basis_set.split("/")[1]
                opt_functional += "-" + opt_dispersion
                opt_functional_and_basis_set = opt_functional + "/" + opt_basis_set
        # Find energies from geometry optimization calculation

        if settings["Boltz energy type"] == "Gibbs free energy":
            e = findall(gibbs_free_energies_regex, data)
            for i in e:
                x = float(i)
                energies.append(x)
                opt_freq_energies.append(x)
        if opt_calc_details and settings["Boltz energy type"] == "Electronic energy":
            data_no_EOL = data.replace("\n ", "")
            if gaussian:
                end_sections = findall(r"1\\1\\.*?\\\\@|1\|1\|.*?\|\|@", data_no_EOL, DOTALL)
                for section in end_sections:
                    section = section.replace("\n ", "")
                    section = section.replace("|", "\\")
                    if sub("^\s+", "", opt_calc_details) in section:  # Find blocks with optimization results
                        e = findall(sp_energies_regex, section)
                        for i in e:
                            x = float(i)
                            energies.append(x)
                            opt_freq_energies.append(x)
            if orca == True:
                geom_opt_energy_blocks = findall(
                    r"\*\*\*\*\*\*\*\*\*\*\*\*HURRAY" r"\*\*\*\*\*\*\*\*\*\*\*\*.*?\*\*\* OPTIMIZATION RUN DONE " 
                    r"\*\*\*",
                    data,
                    DOTALL,
                )
                for block in geom_opt_energy_blocks:
                    energy = findall(sp_energies_regex, block)[0]
                    energies.append(float(energy))
                    opt_freq_energies.append(float(energy))
        # Find vibrational frequencies and IR intensities

        freq_section = findall(freq_section_regex, data, DOTALL)
        for text in freq_section:
            conformer_frequency_groups = findall(frequencies_regex, text)
            conformer_ir_intensity_groups = findall(ir_intensities_regex, text)
            if gaussian:
                conformer_frequency_rotatory_strengths_groups = findall(frequency_rotatory_strengths_regex, text)
                conformer_frequency_dipole_strengths_groups = findall(frequency_dipole_strengths_regex, text)
                conformer_frequencies = []
                conformer_ir_intensities = []
                conformer_frequency_rotatory_strengths = []
                conformer_frequency_dipole_strengths = []
                for index in range(len(conformer_frequency_groups)):
                    conformer_frequencies.append(float(conformer_frequency_groups[index][0]))
                    conformer_frequencies.append(float(conformer_frequency_groups[index][1]))
                    conformer_frequencies.append(float(conformer_frequency_groups[index][2]))
                    conformer_ir_intensities.append(float(conformer_ir_intensity_groups[index][0]))
                    conformer_ir_intensities.append(float(conformer_ir_intensity_groups[index][1]))
                    conformer_ir_intensities.append(float(conformer_ir_intensity_groups[index][2]))
                    if conformer_frequency_rotatory_strengths_groups:
                        conformer_frequency_rotatory_strengths.append(
                            float(conformer_frequency_rotatory_strengths_groups[index][0])
                        )
                        conformer_frequency_rotatory_strengths.append(
                            float(conformer_frequency_rotatory_strengths_groups[index][1])
                        )
                        conformer_frequency_rotatory_strengths.append(
                            float(conformer_frequency_rotatory_strengths_groups[index][2])
                        )
                        conformer_frequency_dipole_strengths.append(
                            float(conformer_frequency_dipole_strengths_groups[index][0])
                        )
                        conformer_frequency_dipole_strengths.append(
                            float(conformer_frequency_dipole_strengths_groups[index][1])
                        )
                        conformer_frequency_dipole_strengths.append(
                            float(conformer_frequency_dipole_strengths_groups[index][2])
                        )
                frequencies.append(conformer_frequencies)
                ir_intensities.append(conformer_ir_intensities)
                if conformer_frequency_rotatory_strengths_groups:
                    frequency_rotatory_strengths.append(conformer_frequency_rotatory_strengths)
                    frequency_dipole_strengths.append(conformer_frequency_dipole_strengths)
            if orca:
                conformer_frequencies = []
                for frequency in conformer_frequency_groups:
                    if frequency != "0.00":
                        frequency = float(frequency)
                        conformer_frequencies.append(frequency)
                frequencies.append(conformer_frequencies)
        # Check for imaginary (negative) frequencies

        for freq in frequencies:
            conf_number = frequencies.index(freq) + 1
            min_frequencies.append(min(freq))
            if min(freq) < 0:
                # Record which conformer(s) had imaginary frequencies

                imaginary_freq_confs.append(conf_number)
        imaginary_freq_confs = list(set(imaginary_freq_confs))

        if orca:
            ir_blocks = findall(orca_ir_section_regex, data, DOTALL)  # Orca excludes imag freqs from this
            # section, but that is OK because this program would exclude any imag freq conformers later anyway,
            # if user doesn't abort analysis altogether.

            for block in ir_blocks:
                conformer_ir_intensities = findall(orca_ir_blocks_regex, block)
                ir_intensities.append(conformer_ir_intensities)

            vcd_blocks = findall(orca_vcd_section_regex, data, DOTALL)

            if vcd_blocks:

                conformer_frequency_rotatory_strengths = []
                for block in vcd_blocks:  # Get frequency rotatory strengths
                    conformer_frequency_rotatory_strengths = findall(frequency_rotatory_strengths_regex, block)
                    frequency_rotatory_strengths.append(conformer_frequency_rotatory_strengths)

                conformer_frequency_extinction_coefficients = []
                for block in ir_blocks:  # Get frequency molar extinction coefficients.

                    conformer_frequency_extinction_coefficients = findall(
                        frequency_dipole_strengths_regex,
                        block
                    )
                    # Note that while Gaussian computes IR dipole strengths, ORCA 6 actually computes IR extinction
                    # coefficients, NOT dipole strengths. This difference is accounted for by SpectroIBIS later in
                    # the final .ir.bil and .docx files.

                    frequency_dipole_strengths.append(conformer_frequency_extinction_coefficients)

        # Find TD-DFT calc theory level and solvation.

        tddft_calc_details = findall(tddft_calc_details_regex, data, IGNORECASE | DOTALL)
        if tddft_calc_details:

            # Remove filenames from parsed calc details. Then, if they no longer contain keyword, remove them.
            if orca:
                for i in tddft_calc_details:
                    unwanted_name = findall(r"%base \".*\"", i[0])
                    x = i[0].replace(unwanted_name[0], "")
                    if 'tddft' not in x.lower():
                        tddft_calc_details.remove(i)
            for i in tddft_calc_details:
                all_tddft_calc_details.append(i)
            if gaussian:
                # Replace tddft_calc_details tuple with new tuple with cleaned data

                c = []
                for i, j in enumerate(tddft_calc_details):
                    b = tddft_calc_details[i][0].replace("\n ", "")
                    a = (b, tddft_calc_details[i][1])
                    c.append(a)
                tddft_calc_details = c
            tddft_solvent = findall(calc_details_solvent_regex, tddft_calc_details[0][0], IGNORECASE)
            if not tddft_solvent and orca:
                tddft_solvent = findall(r'SMDSOLVENT "(\w+-*\w*)"', tddft_calc_details[0][0], IGNORECASE)
            if not tddft_solvent:
                tddft_solvent = ""
            elif tddft_solvent:
                tddft_solvent = tddft_solvent[0]
            if orca:
                keyword_tddft_calc_details = []
                for number, section in enumerate(tddft_calc_details):
                    b = sub(r"\n\| *\d+> +[^!]*", r"", tddft_calc_details[number][0])
                    b = sub("!", "", b)
                    b = sub("NOAUTOSTART|PAL\d+|xyzfile", "", b, IGNORECASE)

                    # Clean up ORCA 6's keywords for Martin's DSD functionals
                    b = sub("REVDSD-PBEP86", "revDSD-PBEP86", b, IGNORECASE)
                    b = sub("REVDOD-PBEP86", "revDOD-PBEP86", b, IGNORECASE)
                    b = sub("DSD-PBEP86/2013", "DSD-PBEP86", b, IGNORECASE)
                    b = sub("DSD-BLYP/2013", "DSD-BLYP", b, IGNORECASE)
                    b = sub("D-PBEP86-D4/2021", "D-PBEP86 D4", b, IGNORECASE)
                    b = sub("D-PBEP86/2021", "D-PBEP86", b, IGNORECASE)
                    keyword_tddft_calc_details.append(b)
                tddft_calc_details = [keyword_tddft_calc_details]
            tddft_dispersion = findall(calc_details_dispersion_regex, tddft_calc_details[0][0], IGNORECASE)
            if not tddft_dispersion:
                tddft_dispersion = ""
            elif tddft_dispersion:
                tddft_dispersion = tddft_dispersion[0]
            if orca:
                # Replace tddft_calc_details tuple with new tuple with cleaned data (i.e. no dispersion correction)

                c = []
                for i, j in enumerate(tddft_calc_details):
                    b = tddft_calc_details[i][0].replace(tddft_dispersion, "")
                    a = (b, tddft_calc_details[i][0])
                    c.append(a)
                tddft_calc_details = c
            tddft_functional_and_basis_set = \
            findall(calc_details_functional_and_basis_set_regex, tddft_calc_details[0][0])[0]
            if orca:
                tddft_functional_and_basis_set = sub(r" +", r"/", tddft_functional_and_basis_set)
            tddft_functional_and_basis_set = sub(
                "#\w? *|opt|EmpiricalDispersion=\w*$|scrf=\w*$|Geom=\w*$|Guess=\w*$", "", tddft_functional_and_basis_set
            )
            #  Reformat functional and basis set

            tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("b3lyp", "B3LYP")
            tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("wB97XD", "ωB97X-D")
            tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("wb97xd", "ωB97X-D")
            tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("wB97MV", "ωB97M-V")
            tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("wb97mv", "ωB97M-V")
            tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("wB", "ωB")
            tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("M062X", "M06-2X")
            tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("m062x", "M06-2X")
            tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("Def2", "def2")
            if "def2" in tddft_functional_and_basis_set and "def2-" not in tddft_functional_and_basis_set:
                tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("def2", "def2-")
            tddft_functional_and_basis_set = tddft_functional_and_basis_set.replace("AUG-cc-", "aug-cc-")
            if tddft_dispersion.casefold() == "gd2":
                tddft_dispersion = "D2"
            elif tddft_dispersion.casefold() == "gd3":
                tddft_dispersion = "D3"
            elif tddft_dispersion.casefold() == "gd3bj":
                tddft_dispersion = "D3(BJ)"
            if tddft_dispersion:
                tddft_functional = tddft_functional_and_basis_set.split("/")[0]
                tddft_basis_set = tddft_functional_and_basis_set.split("/")[1]
                tddft_functional += "-" + tddft_dispersion
                tddft_functional_and_basis_set = tddft_functional + "/" + tddft_basis_set
        # Find NMR calc theory level and solvation.

        nmr_calc_details = findall(nmr_calc_details_regex, data, IGNORECASE | DOTALL)
        if nmr_calc_details:

            # Remove filenames from parsed calc details. Then, if they no longer contain keyword, remove them.
            if orca:
                for i in nmr_calc_details:
                    unwanted_name = findall(r"%base \".*\"", i[0])
                    x = i[0].replace(unwanted_name[0], "")
                    if 'nmr' not in x.lower():
                        nmr_calc_details.remove(i)
            for i in nmr_calc_details:
                all_nmr_calc_details.append(i)
            if gaussian:
                # Replace nmr_calc_details tuple with new tuple with cleaned data

                c = []
                for i, j in enumerate(nmr_calc_details):
                    b = nmr_calc_details[i][0].replace("\n ", "")
                    a = (b, nmr_calc_details[i][1])
                    c.append(a)
                nmr_calc_details = c
            nmr_solvent = findall(calc_details_solvent_regex, nmr_calc_details[0][0], IGNORECASE)
            if not nmr_solvent and orca:
                nmr_solvent = findall(r'SMDSOLVENT "(\w+-*\w*)"', nmr_calc_details[0][0], IGNORECASE)
            if not nmr_solvent:
                nmr_solvent = ""
            elif nmr_solvent:
                nmr_solvent = nmr_solvent[0]
            if orca:
                keyword_nmr_calc_details = []
                for number, section in enumerate(nmr_calc_details):
                    b = sub(r"\n\| *\d+> +[^!]*", r"", nmr_calc_details[number][0])
                    b = sub("!", "", b)
                    b = sub("NOAUTOSTART|PAL\d+|xyzfile", "", b, IGNORECASE)

                    # Clean up ORCA 6's keywords for Martin's DSD functionals
                    b = sub("REVDSD-PBEP86", "revDSD-PBEP86", b, IGNORECASE)
                    b = sub("REVDOD-PBEP86", "revDOD-PBEP86", b, IGNORECASE)
                    b = sub("DSD-PBEP86/2013", "DSD-PBEP86", b, IGNORECASE)
                    b = sub("DSD-BLYP/2013", "DSD-BLYP", b, IGNORECASE)
                    b = sub("D-PBEP86-D4/2021", "D-PBEP86 D4", b, IGNORECASE)
                    b = sub("D-PBEP86/2021", "D-PBEP86", b, IGNORECASE)
                    keyword_nmr_calc_details.append(b)
                nmr_calc_details = [keyword_nmr_calc_details]
            nmr_dispersion = findall(calc_details_dispersion_regex, nmr_calc_details[0][0], IGNORECASE)
            if not nmr_dispersion:
                nmr_dispersion = ""
            elif nmr_dispersion:
                nmr_dispersion = nmr_dispersion[0]
            if gaussian:
                nmr_calc_details = sub("\n ", "", nmr_calc_details[0][0])  # Remove "\n " in case text of interest is
                # split across multiple lines
            if orca:
                nmr_calc_details = sub(nmr_dispersion, "", nmr_calc_details[0][0])
            nmr_functional_and_basis_set = findall(calc_details_functional_and_basis_set_regex, nmr_calc_details)[0]
            if orca:
                nmr_functional_and_basis_set = sub(r" +", r"/", nmr_functional_and_basis_set)
            nmr_functional_and_basis_set = sub(
                "#\w? *|opt|EmpiricalDispersion=\w*$|scrf=\w*$|Geom=\w*$|Guess=\w*$", "", nmr_functional_and_basis_set
            )
            #  Reformat functional and basis set

            nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("b3lyp", "B3LYP")
            nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("wB97XD", "ωB97X-D")
            nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("wb97xd", "ωB97X-D")
            nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("wB97MV", "ωB97M-V")
            nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("wb97mv", "ωB97M-V")
            nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("wB", "ωB")
            nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("M062X", "M06-2X")
            nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("m062x", "M06-2X")
            nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("Def2", "def2")
            if "def2" in nmr_functional_and_basis_set and "def2-" not in nmr_functional_and_basis_set:
                nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("def2", "def2-")
            nmr_functional_and_basis_set = nmr_functional_and_basis_set.replace("AUG-cc-", "aug-cc-")
            if nmr_dispersion.casefold() == "gd2":
                nmr_dispersion = "D2"
            elif nmr_dispersion.casefold() == "gd3":
                nmr_dispersion = "D3"
            elif nmr_dispersion.casefold() == "gd3bj":
                nmr_dispersion = "D3(BJ)"
            if nmr_dispersion:
                nmr_functional = nmr_functional_and_basis_set.split("/")[0]
                nmr_basis_set = nmr_functional_and_basis_set.split("/")[1]
                nmr_functional += "-" + nmr_dispersion
                nmr_functional_and_basis_set = nmr_functional + "/" + nmr_basis_set
        # Find OR calc theory level and solvation.

        or_calc_details = []
        if gaussian:
            or_calc_details = findall(or_calc_details_regex, data, IGNORECASE | DOTALL)
            if or_calc_details:
                for i in or_calc_details:
                    all_or_calc_details.append(i)
                # Replace or_calc_details tuple with new tuple with cleaned data

                c = []
                for i, j in enumerate(or_calc_details):
                    b = or_calc_details[i][0].replace("\n ", "")
                    a = (b, or_calc_details[i][1])
                    c.append(a)
                or_calc_details = c
                or_solvent = findall(calc_details_solvent_regex, or_calc_details[0][0], IGNORECASE)
                if not or_solvent and orca:
                    or_solvent = findall(r'SMDSOLVENT "(\w+-*\w*)"', or_calc_details[0][0], IGNORECASE)
                if not or_solvent:
                    or_solvent = ""
                elif or_solvent:
                    or_solvent = or_solvent[0]
                or_dispersion = findall(calc_details_dispersion_regex, or_calc_details[0][0], IGNORECASE)
                if not or_dispersion:
                    or_dispersion = ""
                elif or_dispersion:
                    or_dispersion = or_dispersion[0]
                or_calc_details = sub("\n ", "",
                                      or_calc_details[0][0])  # Remove "\n " in case text of interest is split
                # across multiple lines

                or_functional_and_basis_set = findall(calc_details_functional_and_basis_set_regex, or_calc_details)[0]
                or_functional_and_basis_set = sub(
                    "#\w? *|opt|EmpiricalDispersion=\w*$|scrf=\w*$|Geom=\w*$|Guess=\w*$", "",
                    or_functional_and_basis_set
                )
                #  Reformat functional and basis set

                or_functional_and_basis_set = or_functional_and_basis_set.replace("b3lyp", "B3LYP")
                or_functional_and_basis_set = or_functional_and_basis_set.replace("wB97XD", "ωB97X-D")
                or_functional_and_basis_set = or_functional_and_basis_set.replace("wb97xd", "ωB97X-D")
                or_functional_and_basis_set = or_functional_and_basis_set.replace("wB97MV", "ωB97M-V")
                or_functional_and_basis_set = or_functional_and_basis_set.replace("wb97mv", "ωB97M-V")
                or_functional_and_basis_set = or_functional_and_basis_set.replace("wB", "ωB")
                or_functional_and_basis_set = or_functional_and_basis_set.replace("M062X", "M06-2X")
                or_functional_and_basis_set = or_functional_and_basis_set.replace("m062x", "M06-2X")
                or_functional_and_basis_set = or_functional_and_basis_set.replace("Def2", "def2")
                if "def2" in or_functional_and_basis_set and "def2-" not in or_functional_and_basis_set:
                    or_functional_and_basis_set = or_functional_and_basis_set.replace("def2", "def2-")
                or_functional_and_basis_set = or_functional_and_basis_set.replace("AUG-cc-", "aug-cc-")
                if or_dispersion.casefold() == "gd2":
                    or_dispersion = "D2"
                elif or_dispersion.casefold() == "gd3":
                    or_dispersion = "D3"
                elif or_dispersion.casefold() == "gd3bj":
                    or_dispersion = "D3(BJ)"
                if or_dispersion:
                    or_functional = or_functional_and_basis_set.split("/")[0]
                    or_basis_set = or_functional_and_basis_set.split("/")[1]
                    or_functional += "-" + or_dispersion
                    or_functional_and_basis_set = or_functional + "/" + or_basis_set
        # Find ECD results section for each conformer

        tddft_section_list = findall(ecd_regex, data, DOTALL)
        # Find calculated excited state transition wavelengths and rotatory strengths (lengths) of each conformer

        for text in tddft_section_list:
            conformer_wavelength_list = findall(wavelength_regex, text)
            conformer_rotatory_strength_list = findall(rotatory_strength_regex, text)
            wavelength_list.append(conformer_wavelength_list)
            rotatory_strength_list.append(conformer_rotatory_strength_list)
            if gaussian:
                conformer_oscillator_strength_list = findall(oscillator_strength_regex, text)
                oscillator_strength_list.append(conformer_oscillator_strength_list)
        if orca:
            uv_section_list = findall(
                r"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n--------------("
                r".*?)ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE "
                r"MOMENTS\n--------------",
                data,
                DOTALL,
            )
            for text in uv_section_list:
                conformer_oscillator_strength_list = findall(oscillator_strength_regex, text)
                oscillator_strength_list.append(conformer_oscillator_strength_list)
        # Find NMR results section for each conformer

        nmr_text_list = findall(nmr_regex, data, DOTALL)
        if nmr_text_list:
            # Find NMR shielding tensors

            for text in nmr_text_list:
                if "[Alpha]D (static) =" in text:  # VCD calcs contain shielding tensors - this avoids these unwanted
                    # shielding tensors being extracted.

                    continue
                conformer_shielding_tensors = findall(shielding_tensors_regex, text)
                shielding_tensors.append(conformer_shielding_tensors)
        # Find optical rotation results section for each conformer

        optrot_section_list = findall(optrot_section_regex, data, DOTALL)
        if optrot_section_list and gaussian:
            # Find optical rotation wavelengths and rotation strengths of each conformer

            for section in optrot_section_list:
                conformer_optrot_wavelengths = []
                conformer_optrot_strengths = []
                conformer_optrot_data = findall(optrot_regex, section, DOTALL)
                for rotation in conformer_optrot_data:
                    wavelength_angstroms = rotation[0]
                    wavelength_nm = str(float(wavelength_angstroms) / 10)
                    conformer_optrot_wavelengths.append(wavelength_nm)
                    optrot_strength = rotation[1]
                    conformer_optrot_strengths.append(optrot_strength)
                optrot_wavelengths.append(conformer_optrot_wavelengths)
                optrot_strengths.append(conformer_optrot_strengths)
        # Find sets of cartesian coordinates

        data_no_EOL = data.replace("\n ", "")
        coordinates_block_list = findall(coords_block_regex, data_no_EOL, DOTALL)
        for conformer_number, text in enumerate(coordinates_block_list):
            text = text.replace("\n ", "")
            conformer = findall(coords_regex, text)
            for atom in conformer:
                conformer_element_list.append(atom[0])
                conformer_x_cartesian_coords_list.append(atom[1])
                conformer_y_cartesian_coords_list.append(atom[2])
                conformer_z_cartesian_coords_list.append(atom[3])
            element_list.append(conformer_element_list)
            x_cartesian_coords_list.append(conformer_x_cartesian_coords_list)
            y_cartesian_coords_list.append(conformer_y_cartesian_coords_list)
            z_cartesian_coords_list.append(conformer_z_cartesian_coords_list)
            conformer_element_list = []
            conformer_x_cartesian_coords_list = []
            conformer_y_cartesian_coords_list = []
            conformer_z_cartesian_coords_list = []
            # Create a log of all coordinates texts for later analysis/checks.

            all_coordinates_texts.append(coordinates_block_list[conformer_number])
        # If required, extract conformer suffixes from .chk filenames in Gaussian output files.

        if gaussian and settings["Mode"] == "Create input files" and conformer_suffix is None:
            a = findall(r"%chk=.*?(conf-\d+|conformer-\d+|M\d\d\d\d)\.chk", data, IGNORECASE)
            for suffix in a:
                if suffix not in chk_conf_suffixes:
                    chk_conf_suffixes.append(suffix)
        # Record which (conjoined) files contained what data in a string, then record this string in a list.

        file_contents = ""
        if coordinates_block_list:
            file_contents += "opt, "
        else:
            file_contents += "no opt, "
        if freq_section:
            file_contents += "freq, "
        else:
            file_contents += "no freq, "
        if e or sp_calc_details:
            file_contents += "energy, "
        else:
            file_contents += "no energy, "
        if nmr_text_list:
            file_contents += "nmr, "
        else:
            file_contents += "no nmr, "
        if tddft_section_list:
            file_contents += "tddft, "
        else:
            file_contents += "no tddft, "
        if or_calc_details:
            file_contents += "optical rotation."
        else:
            file_contents += "no optical rotation."
        list_of_file_contents.append(file_contents)
    # If dedicated single-point energy calculations were detected, use these energies instead of the geom opt energies.

    if sp_functional_and_basis_set:
        energies = []
        if settings["Boltz energy type"] == "Gibbs free energy":
            if len(sp_energies) != len(Gibbs_corrections):
                parser_error_check = "Error detected"
                error_message = (
                        "Extracted "
                        + str(len(sp_energies))
                        + " single-point energies and "
                        + str(len(Gibbs_corrections))
                        + " thermal corrections."
                )
                return (
                    "",
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    parser_error_check,
                    0,
                    error_message,
                )
            for conformer_number, Gibbs_correction in enumerate(Gibbs_corrections):
                energy = sp_energies[conformer_number] + float(Gibbs_correction)
                energies.append(energy)
        elif settings["Boltz energy type"] == "Electronic energy":
            energies = sp_energies
    # Create a filename for the final, analysed data to be saved as files under.

    results_directory = list_of_filepaths[0]
    results_directory = sub(
        "conf-\d+\.log|conformer-\d+\.log|\.log|M\d\d\d\d\.log|conf-\d+\.out|conformer-\d+\.out|M\d\d\d\d\.out|\.out",
        "",
        results_directory,
        IGNORECASE
    )
    results_directory = sub("-$|_$| $", "", results_directory)

    # Perform error checks

    if imaginary_freq_confs:
        for conf_number in imaginary_freq_confs:
            if list_of_conformer_suffixes:
                imag_conf_name = list_of_conformer_suffixes[int(conf_number) - 1]
                if imaginary_freq_confs_text:
                    imaginary_freq_confs_text += ", "
                imaginary_freq_confs_text += imag_conf_name
            else:
                if conf_number == 1:
                    conf_number = str(conf_number) + "st conformer"
                elif conf_number == 2:
                    conf_number = str(conf_number) + "nd conformer"
                elif conf_number == 3:
                    conf_number = str(conf_number) + "rd conformer"
                elif conf_number > 3:
                    conf_number = str(conf_number) + "th conformer"
                if imaginary_freq_confs_text:
                    imaginary_freq_confs_text += ", "
                imaginary_freq_confs_text += str(conf_number)
            parser_error_check = "Imaginary frequency/ies detected"
    if gaussian and settings["Mode"] == "Create input files" and conformer_suffix is None:
        if 0 < len(chk_conf_suffixes) < len(x_cartesian_coords_list):
            parser_error_check = "Error detected"
            error_message = (
                    "Less unique conformer suffixes ("
                    + str(len(chk_conf_suffixes))
                    + ") were extracted from .chk file\nreferences than geometries ("
                    + str(len(x_cartesian_coords_list))
                    + "). Check each conformer has a unique .chk."
            )
        elif len(chk_conf_suffixes) > len(x_cartesian_coords_list):
            parser_error_check = "Error detected"
            error_message = (
                    "More unique conformer suffixes ("
                    + str(len(chk_conf_suffixes))
                    + ") were extracted from .chk file\nreferences than geometries ("
                    + str(len(x_cartesian_coords_list))
                    + ")."
            )
    if len(set(list_of_file_contents)) > 1 and list_of_conformer_suffixes:
        parser_error_check = "Error detected"
        error_message = (
            "Could not match conformer data for all selected files.\nCheck conformer numbers in filenames "
            "match, and all calcs are complete. "
        )
    if len(set(list_of_file_contents)) > 1 and len(list_of_file_contents) != len(set(list_of_file_contents)):
        # Checks if files matched/combined by SpectroIBIS contain different DFT calc types

        if not list_of_conformer_suffixes:
            parser_error_check = "Error detected"
            error_message = (
                "Selected filenames should end with \"conformer-𝘕𝘜𝘔𝘉𝘌𝘙.out\" etc.\nOtherwise, "
                "data can not be matched to the correct conformers. "
            )
        else:
            parser_error_check = "Error detected"
            error_message = (
                "Could not match conformer data for all selected files.\nCheck conformer numbers in "
                "filenames match, and all calcs are complete. "
            )
    # Check all calc details match the number of confs with those data

    if len(x_cartesian_coords_list) < len(all_opt_calc_details):
        parser_error_check = "Error detected"
        error_message = (
            "Extracted too many geom opt theory levels. Check calcs finished\n"
            "and input file(s) don\'t contain \"opt\" letters outside of opt calc section(s)."

        )
    if sp_energies:
        if len(sp_energies) < len(all_sp_calc_details):
            parser_error_check = "Error detected"
            error_message = (
                "Extracted too many SP energy theory levels. Check calcs finished\n"
                "and input file(s) don\'t contain \" sp \" or \" energy \" outside energy calc section(s)."
            )
    if shielding_tensors:
        if len(shielding_tensors) < len(all_nmr_calc_details):
            parser_error_check = "Error detected"
            error_message = (
                "Extracted too many NMR calculation theory levels. Check calcs finished\n"
                "and input file(s) don\'t contain \"nmr\" letters outside of NMR calc section(s)."
            )
    if wavelength_list:
        if len(wavelength_list) < len(all_tddft_calc_details):
            parser_error_check = "Error detected"
            error_message = (
                "Extracted too many TDDFT calculation theory levels. Check calcs finished\n"
                "and input\nfile(s) don\'t contain \"td\(nstates\" or \"tddft\" outside of TDDFT section(s)."
            )
    if optrot_wavelengths:
        if len(optrot_wavelengths) < len(all_or_calc_details):
            parser_error_check = "Error detected"
            error_message = (
                "Extracted too many OR calculation theory levels. Check calcs finished\n"
                "and input\nfile(s) don\'t contain \"Polar=OptRot\" outside of OR calc section(s)."
            )
    if frequencies:
        if max(min_frequencies) < 0:
            parser_error_check = "Error detected"
            error_message = "All conformers have one or more imaginary frequencies!"
    if (
            not len(x_cartesian_coords_list)
            == len(y_cartesian_coords_list)
            == len(z_cartesian_coords_list)
            == len(element_list)
            == len(energies)
            == len(frequencies)
    ):
        parser_error_check = "Error detected"
        error_message = (
            "Unequal amount of opt vs. freq data was extracted for conformers\n"
            "in file(s). Try checking all calcs have finished successfully."
        )
    if (
            not len(x_cartesian_coords_list)
            == len(y_cartesian_coords_list)
            == len(z_cartesian_coords_list)
            == len(element_list)
            == len(energies)
    ):
        parser_error_check = "Error detected"
        error_message = (
            "Unequal amount of energies vs. coordinates data was extracted for\n"
            "conformers in file(s). Try checking all calcs have finished successfully."
        )
    if sp_functional_and_basis_set and energies == opt_freq_energies:
        parser_error_check = "Error detected"
        error_message = (
            "Dedicated single-point energy recalculations were detected,"
            "\nbut could not be incorporated into data analysis."
        )
    if shielding_tensors:
        if (
                len(x_cartesian_coords_list)
                == len(y_cartesian_coords_list)
                == len(z_cartesian_coords_list)
                == len(element_list)
                == len(energies)
                == len(frequencies)
                > len(shielding_tensors)
        ):
            parser_error_check = "Error detected"
            error_message = (
                    "NMR data was extracted for only "
                    + str(len(shielding_tensors))
                    + " / "
                    + str(len(energies))
                    + " conformers. Check\nthat all .out/.log files were selected, and all calculations are complete. "
            )
        elif (
                len(x_cartesian_coords_list)
                == len(y_cartesian_coords_list)
                == len(z_cartesian_coords_list)
                == len(element_list)
                == len(energies)
                == len(frequencies)
                < len(shielding_tensors)
        ):
            parser_error_check = "Error detected"
            error_message = (
                    "Opt freq data was extracted for only "
                    + str(len(energies))
                    + " / "
                    + str(len(shielding_tensors))
                    + " conformers. Check\nthat all .out/.log files were selected, and all calculations are complete. "
            )
        if (
                len(element_list[0])
                > len(shielding_tensors[0])
        ):
            if calc_software == "orca" and findall(
                    "nuclei = ", data, IGNORECASE
            ):  # Checks for ORCA 6 input file keyword to request NMR calculation for only a subset of nuclei
                parser_error_check = "Error detected"
                error_message = (
                        "NMR shielding tensors were extracted for only "
                        + str(len(shielding_tensors[0]))
                        + " / "
                        + str(len(element_list[0]))
                        + " nuclei.\nNMR calculations for subsets of nuclei are not yet supported by SpectroIBIS.\n"
                          "Please recalculate shielding tensors for all nuclei of every conformer."
                )
            else:
                parser_error_check = "Error detected"
                error_message = (
                        "NMR shielding tensors were extracted for only "
                        + str(len(shielding_tensors[0]))
                        + " / "
                        + str(len(element_list[0]))
                        + " nuclei.\nCheck input file(s) contain NMR shielding tensors for all nuclei."
                )
    if wavelength_list:
        if (
                len(x_cartesian_coords_list)
                == len(y_cartesian_coords_list)
                == len(z_cartesian_coords_list)
                == len(element_list)
                == len(energies)
                == len(frequencies)
                > len(wavelength_list)
                == len(rotatory_strength_list)
                == len(oscillator_strength_list)
        ):
            parser_error_check = "Error detected"
            error_message = (
                    "ECD data was extracted for only "
                    + str(len(wavelength_list))
                    + " / "
                    + str(len(energies))
                    + " conformers. Check\nthat all .out/.log files were selected, and all calculations are complete. "
            )
        elif (
                len(x_cartesian_coords_list)
                == len(y_cartesian_coords_list)
                == len(z_cartesian_coords_list)
                == len(element_list)
                == len(energies)
                == len(frequencies)
                < len(wavelength_list)
                == len(rotatory_strength_list)
                == len(oscillator_strength_list)
        ):
            parser_error_check = "Error detected"
            error_message = (
                    "Opt freq data was extracted for only "
                    + str(len(energies))
                    + " / "
                    + str(len(wavelength_list))
                    + " conformers. Check\nthat all .out/.log files were selected, and all calculations are complete. "
            )
        if not len(wavelength_list) == len(rotatory_strength_list) == len(oscillator_strength_list):
            parser_error_check = "Error detected"
            error_message = "Not all required wavelength/ECD/UV data was extracted from file(s)."
    if frequency_rotatory_strengths:
        if (
                len(x_cartesian_coords_list)
                == len(y_cartesian_coords_list)
                == len(z_cartesian_coords_list)
                == len(element_list)
                == len(energies)
                > len(frequencies)
                == len(frequency_rotatory_strengths)
                == len(frequency_dipole_strengths)
        ):
            parser_error_check = "Error detected"
            error_message = (
                    "VCD data was extracted for only "
                    + str(len(frequency_rotatory_strengths))
                    + " / "
                    + str(len(energies))
                    + " conformers. Check\nthat all .out/.log files were selected, and all calculations are complete. "
            )
        elif (
                len(x_cartesian_coords_list)
                == len(y_cartesian_coords_list)
                == len(z_cartesian_coords_list)
                == len(element_list)
                == len(energies)
                < len(frequencies)
                == len(frequency_rotatory_strengths)
                == len(frequency_dipole_strengths)
        ):
            parser_error_check = "Error detected"
            error_message = (
                    "Opt data was extracted for only "
                    + str(len(energies))
                    + " / "
                    + str(len(wavelength_list))
                    + " conformers. Check\nthat all .out/.log files were selected, and all calculations are complete. "
            )
        if not len(frequencies) == len(frequency_rotatory_strengths) == len(frequency_dipole_strengths):
            parser_error_check = "Error detected"
            error_message = "Not all required frequency/VCD/IR data was extracted from file(s)."
    if optrot_wavelengths:
        if (
                len(x_cartesian_coords_list)
                == len(y_cartesian_coords_list)
                == len(z_cartesian_coords_list)
                == len(element_list)
                == len(energies)
                == len(frequencies)
                > len(optrot_wavelengths)
                == len(optrot_strengths)
        ):
            parser_error_check = "Error detected"
            error_message = (
                    "OR data was extracted for only "
                    + str(len(optrot_wavelengths))
                    + " / "
                    + str(len(energies))
                    + " conformers. Check\nthat all .out/.log files were selected, and all calculations are complete. "
            )
        elif (
                len(x_cartesian_coords_list)
                == len(y_cartesian_coords_list)
                == len(z_cartesian_coords_list)
                == len(element_list)
                == len(energies)
                == len(frequencies)
                < len(optrot_wavelengths)
                == len(optrot_strengths)
        ):
            parser_error_check = "Error detected"
            error_message = (
                    "Opt freq data was extracted for only "
                    + str(len(energies))
                    + " / "
                    + str(len(wavelength_list))
                    + " conformers. Check\nthat all .out/.log files were selected, and all calculations are complete. "
            )
        if not len(wavelength_list) == len(rotatory_strength_list) == len(oscillator_strength_list):
            parser_error_check = "Error detected"
            error_message = "Not all required OR wavelength/intensity data was extracted from file(s)."
    if not energies and settings["Boltz energy type"] == "Gibbs free energy":
        parser_error_check = "Error detected"
        error_message = "No Gibbs free energies found in file(s)."
    if not frequencies:
        parser_error_check = "Error detected"
        error_message = (
            "No vibrational frequencies found in file(s).\nPlease include freq calculations for all conformers. "
        )
    if not x_cartesian_coords_list and not y_cartesian_coords_list and not z_cartesian_coords_list:
        parser_error_check = "Error detected"
        error_message = (
            "No geometries found in file(s).\nPlease include opt calculations for all conformers. "
        )
    if not frequencies and not x_cartesian_coords_list and not y_cartesian_coords_list and not z_cartesian_coords_list:
        parser_error_check = "Error detected"
        error_message = (
            "No geometries or vibrational frequencies found in file(s).\nPlease include opt & freq calculations for all conformers. "
        )
    if not energies and settings["Boltz energy type"] == "Electronic energy":
        parser_error_check = "Error detected"
        error_message = "No electronic energies found in file(s)."
    if list_of_conformer_suffixes and are_there_files_without_conf_suffix is True:
        parser_error_check = "Error detected"
        error_message = (
            "Some of the selected filenames have conformer number suffixes\n(e.g., ...conf-12.out), but others do not."
        )
    if len(set(tuple(conformer) for conformer in element_list)) > 1:
        parser_error_check = "Error detected"
        error_message = (
            "Different molecules in selected files, or conformers with\ndifferent atom orders (e.g., C1, C2, H3, O4) in selected files. "
        )
    if (
            not x_cartesian_coords_list
            and not y_cartesian_coords_list
            and not z_cartesian_coords_list
            and not element_list
            and not energies
            and not frequencies
    ):
        parser_error_check = "Error detected"
        error_message = "No opt freq data found in file(s)."
    if (
            not x_cartesian_coords_list
            and not y_cartesian_coords_list
            and not z_cartesian_coords_list
            and not element_list
            and not energies
            and frequencies
    ):
        parser_error_check = "Error detected"
        error_message = "No opt data found in file(s)."
    if (
            shielding_tensors
            and settings["Energies and coordinates table"] is False
            and settings["NMR/ECD/VCD/OR table"] is False
            and settings["NMR csv file"] is False
    ):
        parser_error_check = "Error detected"
        error_message = "No relevant output files are selected in the output settings menu."
    if (
            wavelength_list
            and settings["Energies and coordinates table"] is False
            and settings["NMR/ECD/VCD/OR table"] is False
            and settings["SpecDis .cd.bil file"] is False
            and settings["SpecDis .uv.bil file"] is False
    ):
        parser_error_check = "Error detected"
        error_message = "No relevant output files are selected in the output settings menu."
    if not settings["Temperature (K)"].replace(".", "", 1).isdigit():
        parser_error_check = "Error detected"
        error_message = "The temperature setting contains a non-number value."
    if (
            not settings["Energy cutoff (kcal/mol)"].replace(".", "", 1).isdigit()
            or not settings["MAD cutoff (A)"].replace(".", "", 1).isdigit()
    ):
        parser_error_check = "Error detected"
        error_message = "The redundant conformer cutoff settings contain a non-number value."
    if (
            not settings["H slope"].replace(".", "", 1).replace("-", "", 1).isdigit()
            or not settings["H intercept"].replace(".", "", 1).replace("-", "", 1).isdigit()
            or not settings["C slope"].replace(".", "", 1).replace("-", "", 1).isdigit()
            or not settings["C intercept"].replace(".", "", 1).replace("-", "", 1).isdigit()
    ):
        if (
                not settings["H slope"] == ""
                and settings["H intercept"] == ""
                and settings["C slope"] == ""
                and settings["C intercept"] == ""
        ):
            parser_error_check = "Error detected"
            error_message = "The NMR scaling factor settings contain a non-number value."
    return (
        energies,
        element_list,
        x_cartesian_coords_list,
        y_cartesian_coords_list,
        z_cartesian_coords_list,
        wavelength_list,
        rotatory_strength_list,
        shielding_tensors,
        frequencies,
        coordinates_block_list,
        sp_functional_and_basis_set,
        sp_solvent,
        sp_dispersion,
        opt_functional_and_basis_set,
        opt_solvent,
        opt_dispersion,
        nmr_functional_and_basis_set,
        nmr_solvent,
        nmr_dispersion,
        tddft_functional_and_basis_set,
        tddft_solvent,
        tddft_dispersion,
        all_coordinates_texts,
        results_directory,
        oscillator_strength_list,
        list_of_conformer_suffixes,
        parser_error_check,
        imaginary_freq_confs_text,
        error_message,
        frequency_rotatory_strengths,
        optrot_wavelengths,
        optrot_strengths,
        frequency_dipole_strengths,
        or_functional_and_basis_set,
        or_solvent,
        or_dispersion,
        ir_intensities,
        chk_conf_suffixes,
        list_of_file_contents,
        calc_software,
    )


# Parses XYZ and SDF files, extracts conformer Cartesian coordinates (+ other data for input file creation)
# and checks for redundant conformers and some common (user) errors


def xyz_sdf_parser(list_of_filepaths, settings):
    """Extracts atom coordinates from XYZ/SDF files (as a list of file paths),
    then performs some error checks and removes redundant conformers."""
    # Define lists

    elements = []
    x_coords = []
    y_coords = []
    z_coords = []
    parsed_geometry_data = []
    total_number_atoms = []

    # Set default value for error status

    parser_error_check = "No error detected"
    error_message = ""

    # Create a filename for the final input file(s).

    results_directory = list_of_filepaths[0]
    results_directory = sub("\.xyz|\.sdf", "", results_directory)
    results_directory = sub(
        "conf-\d+\.log|conformer-\d+\.log|\.log|M\d\d\d\d\.log|conf-\d+\.out|conformer-\d+\.out|M\d\d\d\d\.out|\.out",
        "",
        results_directory,
        IGNORECASE
    )

    # Open xyz/sdf file to extract data

    text = ""
    try:
        for filename in list_of_filepaths:
            # Check if dropped xyz file ends with ' Geometries and Energies.xyz'.

            if filename.endswith(" Geometries and Energies.xyz"):
                parser_error_check = "Error detected"
                error_message = (
                    "Conformer order in .xyz file is different to its original .out/.log file.\nUse"
                    " a different .xyz/.sdf/.out/.log file to avoid mismatching conformer data."
                )
                return [], parser_error_check, error_message
            with open(filename, "r") as f:
                new_text = f.read()
                if filename.endswith(".xyz"):
                    number_atoms = findall(r"^[ \t]*(\d+)[ \t]*\n", new_text, MULTILINE)
                    if len(set(number_atoms)) != 1:
                        parser_error_check = "Error detected"
                        error_message = "Inconsistent number of atoms in conformers in .xyz file."
                        return [], parser_error_check, error_message
                    num_atoms = number_atoms[0]
                    conformer_regex = (
                        r"(?:\n[ \t]*\w+[ \t]+-?\d+\.\d+[ \t]+-?\d+\.\d+[ \t]+-?\d+\.\d+){" + str(num_atoms) + "}"
                    )
                    atoms_regex = r"[ \t]*(\w+)[ \t]+(-?\d+\.\d+)[ \t]+(-?\d+\.\d+)[ \t]+(-?\d+\.\d+)"
                if filename.endswith(".sdf"):
                    # Check number of atoms is consistent across conformers

                    number_atoms = findall(r"^[ \t]*(\d+)(?:[ \t]+\d+)+[ \t]+V\d+[ \t]*\n", new_text, MULTILINE)
                    if len(set(number_atoms)) != 1:
                        parser_error_check = "Error detected"
                        error_message = "Inconsistent number of atoms in conformers in .sdf file."
                        return [], parser_error_check, error_message
                    conformer_regex = (
                        r"(?:[ \t]+-?\d+\.\d+[ \t]+-?\d+\.\d+[ \t]+-?\d+\.\d+[ \t]+\w+(?:[ \t]+-?\d+)*\n)+ "
                    )
                    atoms_regex = r"[ \t]+(-?\d+\.\d+)[ \t]+(-?\d+\.\d+)[ \t]+(-?\d+\.\d+)[ \t]+(\w+)"
            for i in number_atoms:
                total_number_atoms.append(i)
            text += new_text
            f.close()
        # Check number of atoms is consistent across all conformers

        if len(set(total_number_atoms)) != 1:
            parser_error_check = "Error detected"
            error_message = (
                "Inconsistent number of atoms in selected .xyz/.sdf file(s).\nCheck selected file(s) "
                "contain the same compound. "
            )
            return [], parser_error_check, error_message
        # Extract conformer geometries

        conformers = findall(conformer_regex, text)
        for index, conformer in enumerate(conformers):
            conf_elements = []
            conf_x_coords = []
            conf_y_coords = []
            conf_z_coords = []
            atoms = findall(atoms_regex, conformer)
            if filename.endswith(".xyz"):
                for atom in atoms:
                    conf_elements.append(atom[0])
                    conf_x_coords.append(atom[1])
                    conf_y_coords.append(atom[2])
                    conf_z_coords.append(atom[3])
            if filename.endswith(".sdf"):
                for atom in atoms:
                    conf_x_coords.append(atom[0])
                    conf_y_coords.append(atom[1])
                    conf_z_coords.append(atom[2])
                    conf_elements.append(atom[3])
            elements.insert(index, conf_elements)
            x_coords.insert(index, conf_x_coords)
            y_coords.insert(index, conf_y_coords)
            z_coords.insert(index, conf_z_coords)
    except:
        return "Error detected", "Unable to read *.xyz / *.sdf file(s).\n Please check for issues in these file(s)."
    a = []
    for h, i in enumerate(elements):
        if h > 0:
            if i != a:  # Check for any conformers with different element lists (e.g. user selected different compounds)
                parser_error_check = "Error detected"
                error_message = (
                    "Inconsistent chemical element lists in selected file(s).\nCheck selected file(s) "
                    "contain the same compound. "
                )
                return [], parser_error_check, error_message
        a = i

    # Define function for building XYZ file format strings for a given conformer number, then creating an RDKit molecule


    def RDKitMoleculeMaker(conformer_number):
        """Function for building XYZ file format strings for a given conformer number,
        then creating an RDKit molecule."""  # For convenience, this code was simply ported from data_analysis.py
        xyz_block = ""
        for i in range(-2, len(elements[conformer_number])):
            if i == -2:
                xyz_block += str(len(elements[conformer_number]))
            elif i == -1:
                xyz_block += "\nMolecule"
            elif i > -1:
                if float(x_coords[conformer_number][i]) >= 0:
                    x = "      " + str(x_coords[conformer_number][i])
                else:
                    x = "     " + str(x_coords[conformer_number][i])
                if float(y_coords[conformer_number][i]) >= 0:
                    y = "      " + str(y_coords[conformer_number][i])
                else:
                    y = "     " + str(y_coords[conformer_number][i])
                if float(z_coords[conformer_number][i]) >= 0:
                    z = "      " + str(z_coords[conformer_number][i])
                else:
                    z = "     " + str(z_coords[conformer_number][i])
                x += "0" * abs(18 - len(x))
                y += "0" * abs(18 - len(y))
                z += "0" * abs(18 - len(z))
                xyz_block += "\n  " + str(elements[conformer_number][i]) + x + y + z
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

    # Check for redundant conformers (purely based on Cartesian coordinates, energies not considered).

    if settings["Check for dup confs in XYZ/SDF files"] is True:
        # Find pairs of conformers with similar or identical geometries

        number_of_confs = len(elements)
        mad_threshold = float(settings["MAD cutoff (A)"])
        duplicate_conformers = []
        dup_conf_numbers = []
        duplicate_conf_details = []
        for conf_number_a, element in enumerate(elements):
            removals = 0
            for conf_number_b in range(number_of_confs):
                conf_number_b = conf_number_b - removals
                if conf_number_b + 1 > len(elements):
                    break
                if conf_number_a < conf_number_b:
                    maximum_atom_deviation = maximum_atom_deviation_calculator(conf_number_a, conf_number_b)
                    if maximum_atom_deviation < mad_threshold:
                        # Record which conformer(s) were duplicates

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
                        if str(duplicated_conf_name).endswith("1") and not str(duplicated_conf_name).endswith("11"):
                            duplicated_conf_name = str(duplicated_conf_name) + "st"
                        elif str(duplicated_conf_name).endswith("2") and not str(duplicated_conf_name).endswith("12"):
                            duplicated_conf_name = str(duplicated_conf_name) + "nd"
                        elif str(duplicated_conf_name).endswith("3") and not str(duplicated_conf_name).endswith("13"):
                            duplicated_conf_name = str(duplicated_conf_name) + "rd"
                        else:
                            duplicated_conf_name = str(duplicated_conf_name) + "th"
                        duplicate_conformers.append(dup_conf_name)
                        if settings["Duplicate conformer details"] is True:
                            dup_conf_text = (
                                str(dup_conf_name)
                                + " conformer is a duplicate of "
                                + str(duplicated_conf_name)
                                + ":\nMAD = "
                                + str(maximum_atom_deviation)
                                + " angstroms\n"
                            )
                            duplicate_conf_details.append(dup_conf_text)
                        # Remove duplicate conformer from data

                        del elements[conf_number_b]
                        del x_coords[conf_number_b]
                        del y_coords[conf_number_b]
                        del z_coords[conf_number_b]
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
        parsed_geometry_data.extend([elements, x_coords, y_coords, z_coords])
        return (
            parsed_geometry_data,
            parser_error_check,
            error_message,
            results_directory,
            duplicate_conf_details,
            duplicate_conformers,
        )
    else:
        parsed_geometry_data.extend([elements, x_coords, y_coords, z_coords])
        return parsed_geometry_data, parser_error_check, error_message, results_directory, [], []
