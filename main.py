# The SpectroIBIS modules are structured as so:
#   PARSERS contains functions which accept a list of filenames (with paths) and return a list of lists containing
#   parsed data by conformer.
#   DATA_ANALYSIS contains a function which accepts this list of lists of parsed data and returns a list of lists of
#   analysed data by conformer.
#   WRITERS accepts this list of lists of analysed data and writes file(s) containing analysed data.
#   MAIN contains the GUI and uses functions imported from the above modules. MAIN also reads/writes user settings
#   from/to a settings text file.


# I am the wisest, for I know that I know nothing - Socrates


# Import modules


import sys
from os import path as os_path, system as os_system
from re import search, findall, sub, DOTALL
from tkinter import (
    Menu,
    Toplevel,
    Label,
    BooleanVar,
    StringVar,
    Checkbutton,
    Radiobutton,
    Entry,
    END,
    Button,
    Text,
    LabelFrame,
    filedialog,
    Canvas,
    Scrollbar,
    Frame,
)
from webbrowser import open_new
from tkinterdnd2 import DND_FILES, TkinterDnD
from data_analysis import analyse
from parsers import parse, xyz_sdf_parser
from writers import (
    nmr_csv_writer,
    ir_csv_writer,
    cd_bil_writer,
    uv_bil_writer,
    vc_bil_writer,
    ir_bil_writer,
    or_bil_writer,
    docx_writer,
    xyz_writer,
    dup_conf_txt_writer,
    input_file_writer,
)


def on_drop(event):
    """Gets dropped filepaths and starts analysis on these selected files."""

    files = event.data
    if files.startswith("{") is True and files.endswith("}") is True:
        files = files.lstrip("{").removesuffix("}")
        list_of_filenames = files.split("} {")
    else:
        list_of_filenames = (
            findall(r"(.+?\.out) ?", files)
            + findall(r"(.+?\.log) ?", files)
            + findall(r"(.+?\.xyz) ?", files)
            + findall(r"(.+?\.sdf) ?", files)
        )
    main(list_of_filenames)


def imag_freq_error_window(settings, list_of_filenames, error_message, imag_freq_conf_list):
    """Opens a dialogue window, triggered by any imaginary frequencies, where the user can decide to proceed
    or abort analysis."""

    def proceed_with_exclusion():
        global program_error_result
        program_error_result = "proceed"
        imag_freq_error_window.destroy()
        return program_error_result

    def abort_run():
        global program_error_result
        program_error_result = "abort"
        imag_freq_error_window.destroy()
        return program_error_result

    def proceed_without_exclusion():
        global program_error_result
        program_error_result = "proceed without exclusion"
        imag_freq_error_window.destroy()
        return program_error_result

    imag_freq_error_window = Toplevel(root)
    imag_freq_error_window.geometry("450x150")
    imag_freq_error_window.title("Error: " + error_message)
    imag_freq_error_window.iconbitmap(icon_path)
    label1 = Label(
        imag_freq_error_window,
        text="One or more imaginary frequencies were detected for:\n" + imag_freq_conf_list + ".",
    )
    label1.config(fg="red")
    label1.pack()
    label2 = Label(imag_freq_error_window, text=" ")
    label2.pack()
    # Look for conformer number suffixes

    conf_suffixes = []
    for filename in list_of_filenames:
        conformer_suffix = search(
            "conf-\d+\.log|conformer-\d+\.log|M\d\d\d\d\.log|conf-\d+\.out|conformer-\d+\.out|M\d\d\d\d\.out",
            filename
        )
        if conformer_suffix != None:
            conformer_suffix = conformer_suffix[0]
            conf_suffixes.append(conformer_suffix)
    if settings["Mode"] == "Create input files" and not conf_suffixes:
        imag_freq_error_window.geometry("450x300")
        Button1 = Button(imag_freq_error_window, text=" Don't proceed (recommended) ", command=abort_run)
        Button2 = Button(imag_freq_error_window, text=" Exclude these conformer(s) ", command=proceed_with_exclusion)
        Button3 = Button(
            imag_freq_error_window,
            text=" Keep these conformer(s) and don't remove redundant or high-energy conformers ",
            command=proceed_without_exclusion,
        )
        Button1.pack()
        Button2.pack()
        Button3.pack()
        label3 = Label(
            imag_freq_error_window,
            text="\nNOTE:\n\n-Excluding non-converged conformers will make any new "
            "input files incompatible\n for combined analysis later with your "
            "current .out/.log file.\n\n-Keeping non-converged conformers will "
            "maintain analysis compatibility.\nAlso, any non-converged "
            "conformers will be excluded later \n(along with redundant conformers)"
            " by SpectroIBIS during the final data analysis.",
            fg="blue",
        )
        label3.pack()
        imag_freq_error_window.protocol("WM_DELETE_WINDOW", abort_run)
        imag_freq_error_window.wait_window(imag_freq_error_window)
    else:
        Button1 = Button(imag_freq_error_window, text=" Don't proceed (recommended) ", command=abort_run)
        Button2 = Button(imag_freq_error_window, text=" Exclude this conformer ", command=proceed_with_exclusion)
        Button1.pack()
        Button2.pack()
        imag_freq_error_window.protocol("WM_DELETE_WINDOW", abort_run)
        imag_freq_error_window.wait_window(imag_freq_error_window)
    return program_error_result


def main(list_of_filenames):
    """Main function of SpectroIBIS - for analysing user-selected files and deciding what to do with them."""

    assert list_of_filenames, "No file(s) selected in file selection window."
    first_file_name = search("/([^/]+)$", list_of_filenames[0]).group(1)
    status_text2 = ""
    status_bar2.config(text=status_text2)
    root.update()
    other_files_text = ""
    if len(list_of_filenames) > 1:
        other_files_text = " and " + str((len(list_of_filenames) - 1)) + " other files"
    # Detect if all files are .xyz or .sdf, or a combination of both.

    xyz_present = 0
    sdf_present = 0
    for filename in list_of_filenames:
        if filename.endswith(".xyz"):
            xyz_present += 1
        if filename.endswith(".sdf"):
            sdf_present += 1
    if xyz_present + sdf_present == len(list_of_filenames) and xyz_present > 0 and sdf_present > 0:  # Mix of xyz and
        # sdf files (bad)

        status_text = (
            "ERROR: Selected files are a mix of .xyz and .sdf files.\nPlease only select one file type at once. "
        )
        status_bar.config(text=status_text, foreground="red")
        return
    # If user has submitted .xyz or .sdf file(s), produce input file(s) instead

    if xyz_present == len(list_of_filenames) or sdf_present == len(list_of_filenames):  # all files are xyz files or
        # all files are sdf files

        suggested_filename = (
            first_file_name.removesuffix(".xyz").removesuffix(".sdf").removesuffix("Geometries and Energies.xyz")
        )
        # Create suggested filename

        suggested_filename = sub(r"-?_?conf-\d+|-?_?conformer-\d+|-?_?M\d\d\d\d", r"", suggested_filename)
        suggested_filename = suggested_filename.removesuffix(" Geometries and Energies")

        if settings["Check for dup confs in XYZ/SDF files"] is True:
            status_text = (
                    "Checking conformers in " + first_file_name + other_files_text + "..."
            )  # bottlenecked by conformer alignment
            status_text2 = "(If this step is really slow, skip it via Settings --> Redundant Conformers)"
            status_bar2.config(text=status_text2, foreground="black")
        else:
            status_text = "Extracting conformers from " + first_file_name + other_files_text + "..."

        status_bar.config(text=status_text, foreground="black")
        root.update()
        parsed_xyz_sdf_data = xyz_sdf_parser(list_of_filenames, settings)
        if parsed_xyz_sdf_data[1] == "Error detected":
            status_text = "ERROR: " + parsed_xyz_sdf_data[2]
            status_bar.config(text=status_text, foreground="red")
            status_bar2.config(text="", foreground="black")
            return
        status_text = "Please enter information for input file(s) in the new window."
        status_bar.config(text=status_text, foreground="black")
        status_bar2.config(text="", foreground="black")
        root.update()
        user_decision = input_file_writer_window(settings, parsed_xyz_sdf_data, suggested_filename, "", False, "", "")
        if user_decision[0] == "Create input file.":  # User has decided to create input file.

            # Write a text file with details about redundant conformers

            duplicate_conf_details = parsed_xyz_sdf_data[4]
            if duplicate_conf_details and settings["Duplicate conformer details"] is True:
                status_text = "Writing redundant conformer details to .txt file..."
                status_bar.config(text=status_text)
                root.update()
                dup_conf_file = dup_conf_txt_writer(parsed_xyz_sdf_data, settings)
                if dup_conf_file[0] == "Error detected":
                    status_text = "ERROR: " + dup_conf_file[1]
                    status_bar.config(text=status_text, foreground="red")
                    return
            duplicate_conformers = parsed_xyz_sdf_data[5]
            if duplicate_conformers:
                if (
                    duplicate_conformers[1].count(",") > 3
                ):  # Avoid overfiling the status bar with many conformer numbers
                    status_text2 = (
                        "Excluded " + duplicate_conformers[0] + " redundant conformer" + duplicate_conformers[2] + "."
                    )
                else:
                    status_text2 = (
                        "Excluded "
                        + duplicate_conformers[0]
                        + " redundant conformer"
                        + duplicate_conformers[2]
                        + " ("
                        + duplicate_conformers[1]
                        + ")."
                    )
                status_bar2.config(text=status_text2)
            elif len(duplicate_conformers) == 0 and settings["Check for dup confs in XYZ/SDF files"] is True:
                status_text2 = "No redundant conformers detected."
                status_bar2.config(text=status_text2)
            status_text = "Writing input file(s)..."
            status_bar.config(text=status_text, foreground="black")
            root.update()
            input_filename = user_decision[1]
            input_file = input_file_writer(parsed_xyz_sdf_data, settings, input_filename)
            if input_file[0] == "Error detected":
                status_text = "ERROR: " + input_file[1]
                status_bar.config(text=status_text, foreground="red")
                status_text2 = ""
                status_bar2.config(text=status_text2)
                return
            plural = ""
            if input_file[3] > 1:
                plural = "s"
            status_text = "Created input file" + plural + ". "
            status_bar.config(text=status_text, fg="green")
            root.update()
            return
        else:  # User has exited window and decided not to create input file.
            status_text = ""
            status_bar.config(text=status_text, fg="black")
            root.update()
            return
    # Extract key data from comp chem output files

    settings["Skip excluding duplicate conformers from input files made from output files"] = False  # Set this setting
    # back to default

    status_text = "Extracting data from " + first_file_name + other_files_text + "..."
    status_bar.config(text=status_text, foreground="black")
    root.update()
    parsed_data = parse(list_of_filenames, settings)
    status_text = "Extracted data for " + str(len(parsed_data[0])) + " conformers..."
    status_bar.config(text=status_text)
    root.update()

    if parsed_data[26] == "Error detected":
        status_text = "ERROR: " + parsed_data[28]
        status_bar.config(text=status_text, foreground="red")
        # error_window(parsed_data[34])

        return
    # Record geometries for possible use in later input file creation, if dubious conformers are to be not removed.

    parsed_elements = parsed_data[1].copy()
    parsed_x_coords = parsed_data[2].copy()
    parsed_y_coords = parsed_data[3].copy()
    parsed_z_coords = parsed_data[4].copy()
    parsed_conf_suffixes = parsed_data[25].copy()
    parsed_chk_conf_suffixes = parsed_data[37].copy()
    parsed_file_content_types = parsed_data[38].copy()

    if parsed_data[26] == "Imaginary frequency/ies detected":
        status_text = "Imaginary frequency detected!"
        status_bar.config(text=status_text, foreground="red")
        imag_freq_error_window_response = imag_freq_error_window(
            settings, list_of_filenames, "Imaginary frequencies", parsed_data[27]
        )
        if imag_freq_error_window_response == "abort":
            if settings["Mode"] == "Analyse output files":
                status_text = "Aborted analysis of " + first_file_name + other_files_text + "."
            if settings["Mode"] == "Create input files":
                status_text = "Aborted creating input files from " + first_file_name + other_files_text + "."
            status_bar.config(text=status_text, foreground="black")
            return
        elif (
            imag_freq_error_window_response == "proceed"
            or imag_freq_error_window_response == "proceed without exclusion"
        ):
            status_bar.config(foreground="black")
    # Process this data

    status_text = "Checking conformers from " + first_file_name + other_files_text + "..."  # main bottleneck
    status_bar.config(text=status_text)
    root.update()
    analysed_data = analyse(parsed_data, settings)
    if analysed_data[33] == "Error detected":
        status_text = "ERROR: " + analysed_data[34]
        status_bar.config(text=status_text, foreground="red")
        return
    # If user has selected the input file creation mode, create input files from .out/.log files instead
    # of analysing them

    if settings["Mode"] == "Create input files":
        status_text = "Please enter information for input file(s) in the new window."
        status_bar.config(text=status_text, foreground="black")
        root.update()
        suggested_filename = sub(
            r"-?_?conf-\d+\.log|-?_?conformer-\d+\.log|-?_?M\d\d\d\d\.log|\.log|-?_?conf-\d+\.out|-?_?conformer-\d"
            r"+\.out|-?_?M\d\d\d\d\.out|\.out",
            r"",
            first_file_name,
        )
        # Create suggested filename
        if len(list_of_filenames) > 1:
            multi_out_file_flag = "Multiple output files"
        else:
            multi_out_file_flag = ""
        if parsed_conf_suffixes == parsed_chk_conf_suffixes == []:
            conformer_renumber_flag = True
        else:
            conformer_renumber_flag = False
        keep_imag_freq_flag = ""
        if "imag_freq_error_window_response" in locals():
            if imag_freq_error_window_response == "proceed without exclusion":
                keep_imag_freq_flag = imag_freq_error_window_response
        user_decision = input_file_writer_window(
            settings,
            analysed_data,
            suggested_filename,
            multi_out_file_flag,
            conformer_renumber_flag,
            keep_imag_freq_flag,
            parsed_file_content_types,
        )
        if user_decision[0] == "Create input file.":  # User has decided to create input file.
            if "imag_freq_error_window_response" in locals():
                if imag_freq_error_window_response == "proceed without exclusion":
                    analysed_data = (
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
                        list_of_filenames[0],
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
                        parsed_elements,
                        parsed_x_coords,
                        parsed_y_coords,
                        parsed_z_coords,
                        parsed_conf_suffixes,
                        "",
                        "",
                        "",
                        parsed_chk_conf_suffixes,
                    )
            if settings["Skip excluding duplicate conformers from input files made from output files"] is True:
                analysed_data = (
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
                    list_of_filenames[0],
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
                    parsed_elements,
                    parsed_x_coords,
                    parsed_y_coords,
                    parsed_z_coords,
                    parsed_conf_suffixes,
                    "",
                    "",
                    "",
                    parsed_chk_conf_suffixes,
                )
                save_new_settings()
            status_text = "Writing input file(s)..."
            status_bar.config(text=status_text, foreground="black")
            root.update()

            # Write a text file with details about redundant conformers

            if settings["Skip excluding duplicate conformers from input files made from output files"] is False:
                duplicate_conf_details = analysed_data[44]
                if duplicate_conf_details and settings["Duplicate conformer details"] is True:
                    status_text = "Writing redundant conformer details to .txt file..."
                    status_bar.config(text=status_text)
                    root.update()
                    dup_conf_file = dup_conf_txt_writer(analysed_data, settings)
                    if dup_conf_file[0] == "Error detected":
                        status_text = "ERROR: " + dup_conf_file[1]
                        status_bar.config(text=status_text, foreground="red")
                        return
            input_filename = user_decision[1]
            input_file = input_file_writer(analysed_data, settings, input_filename)
            if input_file[0] == "Error detected":
                status_text = "ERROR: " + input_file[1]
                status_bar.config(text=status_text, foreground="red")
                return
            number_high_energy_confs_removed = input_file[2]
            rel_e_threshold_unit = ""
            if str(settings["Geometry reopt relative energy unit"]) == "kJmol-1":
                rel_e_threshold_unit = "kJ/mol"
            elif str(settings["Geometry reopt relative energy unit"]) == "kcalmol-1":
                rel_e_threshold_unit = "kcal/mol"
            if settings["Geometry reopt"] is True:
                if number_high_energy_confs_removed == 0:
                    removals_text = (
                        "No conformers found above "
                        + str(settings["Geometry reopt relative energy threshold"])
                        + " "
                        + rel_e_threshold_unit
                        + "."
                    )
                elif number_high_energy_confs_removed == 0:
                    removals_text = (
                        "Excluded 1 conformer above "
                        + str(settings["Geometry reopt relative energy threshold"])
                        + " "
                        + rel_e_threshold_unit
                        + "."
                    )
                else:
                    removals_text = (
                        "Excluded "
                        + str(number_high_energy_confs_removed)
                        + " conformers above "
                        + str(settings["Geometry reopt relative energy threshold"])
                        + " "
                        + rel_e_threshold_unit
                        + "."
                    )
            else:
                removals_text = ""
            plural = ""
            if input_file[3] > 1:
                plural = "s"
            status_text = "Created input file" + plural + ". " + removals_text
            status_bar.config(text=status_text, fg="green")
            root.update()

            duplicate_conformers = analysed_data[32]
            if duplicate_conformers:
                if (
                    duplicate_conformers[1].count(",") > 3
                ):  # Avoids over-filling the status bar with many conformer numbers
                    status_text2 = (
                        "Excluded " + duplicate_conformers[0] + " redundant conformer" + duplicate_conformers[2] + "."
                    )
                else:
                    status_text2 = (
                        "Excluded "
                        + duplicate_conformers[0]
                        + " redundant conformer"
                        + duplicate_conformers[2]
                        + " ("
                        + duplicate_conformers[1]
                        + ")."
                    )
                status_bar2.config(text=status_text2)
                root.update()
            elif not duplicate_conformers:
                status_text2 = "No redundant conformers detected." # TODO Below x kj/mol?
                status_bar2.config(text=status_text2)
                root.update()
            return
        else:  # User has exited window and decided not to create input file.
            status_text = "Exited input file template window."
            status_bar.config(text=status_text, fg="black")
            status_text2 = "No input file(s) created."
            status_bar2.config(text=status_text2, fg="black")
            root.update()
            return
    # Make csv file with Boltzmann-averaged shielding tensors

    boltz_shielding_tensors = analysed_data[10]
    if boltz_shielding_tensors and settings["NMR csv file"] is True:
        status_text = "Writing NMR results to csv file..."
        status_bar.config(text=status_text)
        root.update()
        nmr_csv = nmr_csv_writer(analysed_data)
        if nmr_csv[0] == "Error detected":
            status_text = "ERROR: " + nmr_csv[1]
            status_bar.config(text=status_text, foreground="red")
            return
    # Make csv file with Boltzmann-averaged frequencies

    boltz_frequencies = analysed_data[29]
    if boltz_frequencies and settings["Freq csv file"] is True:
        status_text = "Writing frequencies to csv file..."
        status_bar.config(text=status_text)
        root.update()
        freq_csv = ir_csv_writer(analysed_data, settings)
        if freq_csv[0] == "Error detected":
            status_text = "ERROR: " + freq_csv[1]
            status_bar.config(text=status_text, foreground="red")
            return
    # Make .cd.bil file with Boltzmann-averaged ECD data for SpecDis

    boltz_wavelengths = analysed_data[13]
    if boltz_wavelengths and settings["SpecDis .cd.bil file"] is True:
        status_text = "Writing ECD results to SpecDis .cd.bil file..."
        status_bar.config(text=status_text)
        root.update()
        cd_bil = cd_bil_writer(analysed_data)
        if cd_bil[0] == "Error detected":
            status_text = "ERROR: " + cd_bil[1]
            status_bar.config(text=status_text, foreground="red")
            return
    # Make .uv.bil file with Boltzmann-averaged UV data for SpecDis

    boltz_wavelengths = analysed_data[13]
    if boltz_wavelengths and settings["SpecDis .uv.bil file"] is True:
        status_text = "Writing UV results to SpecDis .uv.bil file..."
        status_bar.config(text=status_text)
        root.update()
        uv_bil = uv_bil_writer(analysed_data)
        if uv_bil[0] == "Error detected":
            status_text = "ERROR: " + uv_bil[1]
            status_bar.config(text=status_text, foreground="red")
            return
    # Make .vc.bil file with Boltzmann-averaged VCD data for SpecDis

    boltz_frequency_rotatory_strengths = analysed_data[35]
    if boltz_frequency_rotatory_strengths and settings["SpecDis .vc.bil file"] is True:
        status_text = "Writing VCD results to SpecDis .vc.bil file..."
        status_bar.config(text=status_text)
        root.update()
        vc_bil_writer(analysed_data)
        if vc_bil_writer(analysed_data)[0] == "Error detected":
            status_text = "ERROR: " + uv_bil_writer(analysed_data)[1]
            status_bar.config(text=status_text, foreground="red")
            return
    # Make .ir.bil file with Boltzmann-averaged IR data for SpecDis

    boltz_frequency_dipole_strengths = analysed_data[40]
    if boltz_frequency_dipole_strengths and settings["SpecDis .ir.bil file"] is True:
        status_text = "Writing IR results to SpecDis .ir.bil file..."
        status_bar.config(text=status_text)
        root.update()
        ir_bil = ir_bil_writer(analysed_data)
        if ir_bil[0] == "Error detected":
            status_text = "ERROR: " + ir_bil[1]
            status_bar.config(text=status_text, foreground="red")
            return
    # Make .or.bil file with Boltzmann-averaged VCD data for SpecDis

    ordered_optrot_wavelengths_list = analysed_data[36]
    if ordered_optrot_wavelengths_list and settings["SpecDis .or.bil file"] is True:
        status_text = "Writing optical rotation results to SpecDis .or.bil file..."
        status_bar.config(text=status_text)
        root.update()
        or_bil = or_bil_writer(analysed_data)
        if or_bil[0] == "Error detected":
            status_text = "ERROR: " + or_bil[1]
            status_bar.config(text=status_text, foreground="red")
            return
    # Make a .xyz file containing conformer energies and Cartesian coordinates

    if settings["Write .xyz file"] is True:
        status_text = "Writing conformer geometries and energies to .xyz file..."
        status_bar.config(text=status_text)
        root.update()
        xyz_file = xyz_writer(analysed_data)
        if xyz_file[0] == "Error detected":
            status_text = "ERROR: " + xyz_file[1]
            status_bar.config(text=status_text, foreground="red")
            return
    # Write a text file with details about redundant conformers

    duplicate_conf_details = analysed_data[44]
    if duplicate_conf_details and settings["Duplicate conformer details"] is True:
        status_text = "Writing redundant conformer details to .txt file..."
        status_bar.config(text=status_text)
        root.update()
        dup_conf_file = dup_conf_txt_writer(analysed_data, settings)
        if dup_conf_file[0] == "Error detected":
            status_text = "ERROR: " + dup_conf_file[1]
            status_bar.config(text=status_text, foreground="red")
            return
    # Write a Word document (.docx) with supplementary data

    status_text = "Writing Word document..."
    status_bar.config(text=status_text)
    duplicate_conformers = analysed_data[32]
    if duplicate_conformers:
        if duplicate_conformers[1].count(",") > 3:  # Avoid overfiling the status bar with many conformer numbers
            status_text2 = "Excluded " + duplicate_conformers[0] + " redundant conformer" + duplicate_conformers[
                2] + "."
        else:
            status_text2 = (
                    "Excluded "
                    + duplicate_conformers[0]
                    + " redundant conformer"
                    + duplicate_conformers[2]
                    + " ("
                    + duplicate_conformers[1]
                    + ")."
            )
        status_bar2.config(text=status_text2)
    elif not duplicate_conformers:
        status_text2 = "No redundant conformers detected."
        status_bar2.config(text=status_text2)
    root.update()
    docx_file = docx_writer(analysed_data, settings)
    if docx_file[0] == "Error detected":
        status_text = "ERROR: " + docx_file[1]
        status_bar.config(text=status_text, foreground="red")
        status_text2 = ""
        status_bar2.config(text=status_text2)
        return
    results_name_and_directory = analysed_data[27]
    results_name = search("/([^/]+)$", results_name_and_directory).group(1)
    status_text = "Finished for " + first_file_name + other_files_text + "."
    status_bar.config(text=status_text, fg="green")
    root.update()


def save_new_settings():
    """Saves new settings to the settings text file."""
    if not os_path.isdir(files_folder_path):  # Check required folder exists
        raise FileNotFoundError(
            "SpectroIBIS.exe could not find its associated folder \"Files for SpectroIBIS\". Try moving this "
            "SpectroIBIS.exe file back into its original unzipped folder, "
            "or instead try moving the \"Files for SpectroIBIS\" folder next to this SpectroIBIS.exe file."
        )
    with open(settings_path, "w", encoding="utf-8") as settings_file:
        for key, value in settings.items():
            if key != "Input File Texts":
                settings_file.write(key + ": " + str(value) + "\n")
            if key == "Input File Texts":
                settings_file.write("Input file template is below:")
                for calc_number in range(0, len(settings["Input File Texts"])):
                    if calc_number == 0:
                        settings_file.write("\n*#*#*#* Multi Job Text START *#*#*#*")
                        settings_file.write("\n" + settings["Input File Texts"][calc_number])
                        settings_file.write("\n*#*#*#* Multi Job Text END *#*#*#*")
                    else:
                        settings_file.write("\n*#*#*#* Calculation #" + str(calc_number) + " START *#*#*#*")
                        settings_file.write("\n" + settings["Input File Texts"][calc_number])
                        settings_file.write("\n*#*#*#* Calculation #" + str(calc_number) + " END *#*#*#*")


def read_settings():
    """Looks for and loads settings from the settings text file (if none, makes a new one with default settings)."""

    default_settings = {
        "Mode": "Analyse output files",
        "Energies and coordinates table": True,
        "NMR/ECD/VCD/OR table": True,
        "UV in ECD table": False,
        "NMR csv file": True,
        "Freq csv file": False,
        "SpecDis .cd.bil file": True,
        "SpecDis .uv.bil file": True,
        "SpecDis .vc.bil file": True,
        "SpecDis .ir.bil file": True,
        "SpecDis .or.bil file": True,
        "Write .xyz file": False,
        "Relative energy unit": "kJ/mol",
        "Duplicate conformer details": True,
        "Check for dup confs in XYZ/SDF files": True,
        "Geometry reopt relative energy unit": "kJmol-1",
        "Geometry reopt relative energy threshold": "40.0",
        "Geometry reopt": True,
        "Energy cutoff (kcal/mol)": "0.1",
        "MAD cutoff (A)": "0.1",
        "H slope": "",
        "H intercept": "",
        "C slope": "",
        "C intercept": "",
        "IR freq scaling factor": "",
        "Temperature (K)": "298.15",
        "Boltz energy type": "Gibbs free energy",
        "Input File Conformers Together": False,
        "Skip excluding duplicate conformers from input files made from output files": False,
        "Input File Texts": [
            "\n--Link1--",
            "%chk=⫷⫷⫷COMPOUND NAME⫸⫸⫸_conf-⫷⫷⫷CONFORMER NUMBER⫸⫸⫸.chk "
            "\n%nproc=4\n%mem=50000MB\n#n wB97XD/def2TZVP Opt Freq=NoRaman scrf=("
            "solvent=acetonitrile)\n\n⫷⫷⫷COMPOUND NAME⫸⫸⫸ Conformer "
            "#⫷⫷⫷CONFORMER NUMBER⫸⫸⫸ opt freq calc in acetonitrile\n\n0 1\n⫷⫷⫷ATOM "
            "COORDINATES⫸⫸⫸",
            "%chk=⫷⫷⫷COMPOUND NAME⫸⫸⫸_conf-⫷⫷⫷CONFORMER NUMBER⫸⫸⫸.chk "
            "\n%nproc=4\n%mem=50000MB\n# td(NStates=40) wB97XD/def2TZVP scrf=("
            "solvent=acetonitrile) Geom=Check Guess=Read\n\n⫷⫷⫷COMPOUND NAME⫸⫸⫸ "
            "Conformer #⫷⫷⫷CONFORMER NUMBER⫸⫸⫸ ECD calc in acetonitrile\n\n0 1",
        ],
    }
    global settings
    settings = {}
    input_file_text_mode = False
    if os_path.isfile(settings_path):  # Try to open a pre-existing settings file
        with open(settings_path, "r", encoding="utf-8") as settings_file:
            text = settings_file.read()
            lines = text.split("\n")
            for line in lines:
                if line != "":
                    if line == "Input file template is below:":
                        input_file_text_mode = True
                        break
                    setting = line.split(": ")
                    if setting[1] == "True":
                        setting[1] = True
                    if setting[1] == "False":
                        setting[1] = False
                    settings[setting[0]] = setting[1]
            if input_file_text_mode is True:
                multi_job_text = findall(
                    r"\*#\*#\*#\* Multi Job Text START \*#\*#\*#\*\n(.*?)\n\*#\*#\*#\* Multi Job Text END \*#\*#\*#\*",
                    text,
                    DOTALL,
                )
                input_file_calcs = findall(
                    r"\*#\*#\*#\* Calculation #\d+ START \*#\*#\*#\*\n(.*?)\n\*#\*#\*#\* Calculation #\d+ END "
                    r"\*#\*#\*#\*",
                    text,
                    DOTALL,
                )
                settings["Input File Texts"] = multi_job_text + input_file_calcs
    elif not os_path.isfile(settings_path):  # Write default settings to new settings file, then read it.
        settings = default_settings
        save_new_settings()
        with open(settings_path, "r", encoding="utf-8") as settings_file:
            whole_text = settings_file.read()
            text = settings_file.read().split("\n")
            for line in text:
                if line != "":
                    if line == "Input file template is below:":
                        input_file_text_mode = True
                        break
                    setting = line.split(": ")
                    settings[setting[0]] = setting[1]
            if input_file_text_mode is True:
                multi_job_text = findall(
                    r"\*#\*#\*#\* Multi Job Text START \*#\*#\*#\*\n(.*?)\*#\*#\*#\* Multi Job Text END \*#\*#\*#\*",
                    whole_text,
                    DOTALL,
                )
                input_file_calcs = findall(
                    r"\*#\*#\*#\* Calculation #\d+ START \*#\*#\*#\*\n(.*?)\*#\*#\*#\* Calculation #\d+ END "
                    r"\*#\*#\*#\*",
                    whole_text,
                    DOTALL,
                )
                settings["Input File Texts"] = multi_job_text + input_file_calcs
    return settings


def file_selection():
    """Adds a file opener window to file menu, which can start analysis on user-selected OUT/XYZ/SDF files.
    This is an alternative file selection method to drag-and-drop."""

    global list_of_filenames
    list_of_filenames = filedialog.askopenfilenames(
        title="Select Files",
        filetypes=[("Output File(s)", ".out .log"), ("XYZ File(s)", ".xyz"), ("SDF File(s)", ".sdf")],
    )
    list_of_filenames = list(list_of_filenames)
    main(list_of_filenames)


# Determine if application is running as a script file or frozen executable file.
# This is useful for finding program files when launched via a shortcut.


if getattr(sys, "frozen", False):
    application_path = os_path.dirname(sys.executable)
elif __file__:
    application_path = os_path.dirname(__file__)
# Set file paths for icon, manual and settings.


global icon_path
global manual_path
global settings_path
files_folder_path = os_path.join(application_path, "Files For SpectroIBIS")
icon_path = os_path.join(files_folder_path, "SpectroIBIS_icon.ico")
manual_path = os_path.join(files_folder_path, "SpectroIBIS Manual.pdf")
settings_path = os_path.join(files_folder_path, "SpectroIBIS Settings.txt")
settings = read_settings()
# Launch GUI


root = TkinterDnD.Tk()
root.resizable(False, False)
root.title("SpectroIBIS")
root.geometry("400x363")
try:
    root.iconbitmap(icon_path)
except:
    raise FileNotFoundError(
        "SpectroIBIS.exe could not find its icon file \"SpectroIBIS_icon.ico\" in a neighboring folder \"Files for "
        "SpectroIBIS\". Try moving this SpectroIBIS.exe file back into its original unzipped folder."
    )
# Create a menu bar


menu_bar = Menu(root)
root.config(menu=menu_bar)

# Create a file menu


file_menu = Menu(menu_bar, tearoff=0)
menu_bar.add_cascade(label="File", menu=file_menu)
file_menu.add_command(label="Open", command=file_selection)


# Create a box for user to drag and drop their files to


drop_box = Text(root, height=14, width=48)
drop_box.pack()
if settings["Mode"] == "Analyse output files":
    drop_box_text = """

To analyse data:
Drag and drop Gaussian/ORCA output files for
all conformers of your compound here, together.


OR


To create input files:
Drag and drop XYZ or SDF files instead.

"""
else:
    drop_box_text = """





To create input files:
Drag and drop XYZ, SDF or output files here.





"""
drop_box.insert(END, drop_box_text, "center")
drop_box.tag_configure("center", justify="center")
drop_box.config(state="disabled")
drop_box.configure(font=("Segoe UI", 12))
drop_box.configure()


def open_manual():
    """Opens the user manual PDF file."""

    # Make manual_path have quotations around directories with a space in name

    a = ""
    for dir in manual_path.split("\\"):
        if " " in dir:
            dir = '"' + dir + '"'
        a += "\\" + dir
    a = a.removeprefix("\\")

    os_system(a)  # Open manual


def open_about_window():
    """Opens a window with general information about SpectroIBIS."""

    window = Toplevel(root)
    window.title("About")
    window.iconbitmap(icon_path)
    label_frame = LabelFrame(window, borderwidth=0, highlightthickness=0)
    label_frame.pack()
    label0 = Label(label_frame, text="Spectro", font=("Colonna MT", 36, "italic"))
    label0.grid(row=0, column=0)
    label00 = Label(label_frame, text="IBIS", font=("Colonna MT", 36, "bold"))
    label00.grid(row=0, column=1)
    label1 = Label(window, text="Version 1.0.0", font=("Arial", 14))
    label1.pack()
    label2 = Label(
        window, text="https://github.com/bbulcock/SpectroIBIS", font=("Arial", 10), fg="blue", cursor="hand2"
    )
    label2.pack()
    label2b = Label(
        window, text="Automated data processing for quantum chemical spectroscopic calculations", font=("Arial", 10)
    )
    label2b.pack()
    label2.bind("<Button-1>", lambda e: open_new("https://github.com/bbulcock/SpectroIBIS"))
    label3 = Label(window, text="Developed by Brodie Bulcock", font=("Arial", 10))
    label3.pack()
    label4 = Label(window, text="\nIf you found this tool useful, please cite:\n(Publication coming soon!)",
                   font=("Arial", 10))
    label4.pack()


def output_window():
    """Opens a window where the user can modify settings for SpectroIBIS data output."""

    def get_parameters():
        """Retrieves newly entered settings and saves these to the settings text file."""
        settings["Energies and coordinates table"] = var1.get()
        settings["NMR/ECD/VCD/OR table"] = var4.get()
        settings["UV in ECD table"] = var4b.get()
        settings["NMR csv file"] = var5.get()
        settings["SpecDis .cd.bil file"] = var6.get()
        settings["SpecDis .uv.bil file"] = var7.get()
        settings["SpecDis .vc.bil file"] = var8.get()
        settings["SpecDis .ir.bil file"] = var9.get()
        settings["SpecDis .or.bil file"] = var10.get()
        settings["Relative energy unit"] = var11.get()
        settings["Write .xyz file"] = var12.get()
        settings["Freq csv file"] = var13.get()
        settings["Mode"] = var14.get()
        if var14.get() == "Create input files":
            checkbox1.config(state="disabled")
            checkbox4.config(state="disabled")
            checkbox4b.config(state="disabled")
            checkbox5.config(state="disabled")
            checkbox6.config(state="disabled")
            checkbox6.config(state="disabled")
            checkbox7.config(state="disabled")
            checkbox8.config(state="disabled")
            checkbox9.config(state="disabled")
            checkbox10.config(state="disabled")
            checkbox11.config(state="disabled")
            checkbox12.config(state="disabled")
            radiobutton1.config(state="disabled")
            radiobutton2.config(state="disabled")
            boltz_data_frame.config(fg="grey")
            ms_word_frame.config(fg="grey")
            other_files_frame.config(fg="grey")
        else:
            checkbox1.config(state="normal")
            checkbox4.config(state="normal")
            checkbox4b.config(state="normal")
            checkbox5.config(state="normal")
            checkbox6.config(state="normal")
            checkbox6.config(state="normal")
            checkbox7.config(state="normal")
            checkbox8.config(state="normal")
            checkbox9.config(state="normal")
            checkbox10.config(state="normal")
            checkbox11.config(state="normal")
            checkbox12.config(state="normal")
            radiobutton1.config(state="normal")
            radiobutton2.config(state="normal")
            boltz_data_frame.config(fg="black")
            ms_word_frame.config(fg="black")
            other_files_frame.config(fg="black")
        # Update main drop box text if settings["Mode"] is changed.

        if settings["Mode"] == "Analyse output files":
            drop_box_text = """

To analyse data:
Drag and drop Gaussian/ORCA output files for
all conformers of your compound here, together.


OR


To create input files:
Drag and drop XYZ or SDF files instead.

"""
        else:
            drop_box_text = """





To create input files:
Drag and drop XYZ, SDF or output files here.





"""
        drop_box.config(state="normal")
        drop_box.delete("1.0", END)
        drop_box.insert(END, drop_box_text, "center")
        drop_box.tag_configure("center", justify="center")
        drop_box.config(state="disabled")
        drop_box.update()
        save_new_settings()
        status_text = ""
        status_bar.config(text=status_text)
        status_text2 = ""
        status_bar2.config(text=status_text2)

    window = Toplevel(root)
    window.title("Output Settings")
    window.iconbitmap(icon_path)
    var1 = BooleanVar()
    var4 = BooleanVar()
    var4b = BooleanVar()
    var5 = BooleanVar()
    var6 = BooleanVar()
    var7 = BooleanVar()
    var8 = BooleanVar()
    var9 = BooleanVar()
    var10 = BooleanVar()
    var11 = StringVar()
    var12 = BooleanVar()
    var13 = BooleanVar()
    var14 = StringVar()
    var1.set(settings["Energies and coordinates table"])
    var4.set(settings["NMR/ECD/VCD/OR table"])
    var4b.set(settings["UV in ECD table"])
    var5.set(settings["NMR csv file"])
    var6.set(settings["SpecDis .cd.bil file"])
    var7.set(settings["SpecDis .uv.bil file"])
    var8.set(settings["SpecDis .vc.bil file"])
    var9.set(settings["SpecDis .ir.bil file"])
    var10.set(settings["SpecDis .or.bil file"])
    var11.set(settings["Relative energy unit"])
    var12.set(settings["Write .xyz file"])
    var13.set(settings["Freq csv file"])
    var14.set(settings["Mode"])

    mode_frame = LabelFrame(window, text="What to do with Gaussian/ORCA output files")
    mode_frame.grid(row=0, column=1, sticky="w")
    space0 = Label(window, text="")
    space0.grid(row=1, column=1, sticky="w")
    boltz_data_frame = LabelFrame(window, text="Boltzmann-Averaged Data (When Available)")
    boltz_data_frame.grid(row=2, column=1, sticky="w")
    space1 = Label(window, text="")
    space1.grid(row=3, column=1, sticky="w")
    ms_word_frame = LabelFrame(window, text="Supplementary Data (.docx file)")
    ms_word_frame.grid(row=4, column=1, sticky="w")
    space2 = Label(window, text="")
    space2.grid(row=5, column=1, sticky="w")
    other_files_frame = LabelFrame(window, text="Other")
    other_files_frame.grid(row=6, column=1, sticky="w")
    space3 = Label(window, text="")
    space3.grid(row=7, column=1, sticky="w")

    radiobutton_a = Radiobutton(
        mode_frame,
        text="Analyse data",
        variable=var14,
        value="Analyse output files",
        anchor="w",
        command=get_parameters,
    )
    radiobutton_a.pack(side="left", anchor="center", padx=34)
    radiobutton_b = Radiobutton(
        mode_frame,
        text="Create input file(s) ",
        variable=var14,
        value="Create input files",
        anchor="w",
        command=get_parameters,
    )
    radiobutton_b.pack(side="right", anchor="center", padx=34)
    checkbox1 = Checkbutton(
        ms_word_frame, text="Energies and geometries table", variable=var1, anchor="w", command=get_parameters
    )
    checkbox1.grid(row=0, column=1, sticky="w")
    radiobutton1 = Radiobutton(
        ms_word_frame,
        text="Relative energies in kcal/mol",
        variable=var11,
        value="kcal/mol",
        anchor="w",
        command=get_parameters,
    )
    radiobutton1.grid(row=1, column=1, sticky="w")
    radiobutton2 = Radiobutton(
        ms_word_frame,
        text="Relative energies in kJ/mol  ",
        variable=var11,
        value="kJ/mol",
        anchor="w",
        command=get_parameters,
    )
    radiobutton2.grid(row=1, column=2, sticky="w")
    checkbox4 = Checkbutton(
        ms_word_frame, text="NMR/ECD/VCD/OR table", variable=var4, anchor="w", command=get_parameters
    )
    checkbox4.grid(row=0, column=2, sticky="w")
    checkbox4b = Checkbutton(
        ms_word_frame, text="Include UV in ECD table", variable=var4b, anchor="w", command=get_parameters
    )
    checkbox4b.grid(row=2, column=1, sticky="w")
    checkbox5 = Checkbutton(boltz_data_frame, text="NMR csv file", variable=var5, anchor="w", command=get_parameters)
    checkbox5.grid(row=3, column=1, sticky="w")
    checkbox6 = Checkbutton(
        boltz_data_frame, text="SpecDis cd.bil file", variable=var6, anchor="w", command=get_parameters
    )
    checkbox6.grid(row=3, column=3, sticky="w")
    checkbox7 = Checkbutton(
        boltz_data_frame, text="SpecDis uv.bil file", variable=var7, anchor="w", command=get_parameters
    )
    checkbox7.grid(row=4, column=1, sticky="w")
    checkbox8 = Checkbutton(
        boltz_data_frame, text="SpecDis vc.bil file", variable=var8, anchor="w", command=get_parameters
    )
    checkbox8.grid(row=4, column=2, sticky="w")
    checkbox9 = Checkbutton(
        boltz_data_frame, text="SpecDis ir.bil file", variable=var9, anchor="w", command=get_parameters
    )
    checkbox9.grid(row=4, column=3, sticky="w")
    checkbox10 = Checkbutton(
        boltz_data_frame, text="SpecDis or.bil file", variable=var10, anchor="w", command=get_parameters
    )
    checkbox10.grid(row=5, column=1, sticky="w")
    checkbox11 = Checkbutton(
        other_files_frame,
        text="Geometries and energies in xyz file, ordered by energy",
        variable=var12,
        anchor="w",
        command=get_parameters,
    )
    checkbox11.grid(row=5, column=1, sticky="w")
    checkbox12 = Checkbutton(boltz_data_frame, text="IR csv file", variable=var13, anchor="w", command=get_parameters)
    checkbox12.grid(row=3, column=2, sticky="w")

    # Update relevant widgets

    if var14.get() == "Create input files":
        checkbox1.config(state="disabled")
        checkbox4.config(state="disabled")
        checkbox4b.config(state="disabled")
        checkbox5.config(state="disabled")
        checkbox6.config(state="disabled")
        checkbox6.config(state="disabled")
        checkbox7.config(state="disabled")
        checkbox8.config(state="disabled")
        checkbox9.config(state="disabled")
        checkbox10.config(state="disabled")
        checkbox11.config(state="disabled")
        checkbox12.config(state="disabled")
        radiobutton1.config(state="disabled")
        radiobutton2.config(state="disabled")
        boltz_data_frame.config(fg="grey")
        ms_word_frame.config(fg="grey")
        other_files_frame.config(fg="grey")
    else:
        checkbox1.config(state="normal")
        checkbox4.config(state="normal")
        checkbox4b.config(state="normal")
        checkbox5.config(state="normal")
        checkbox6.config(state="normal")
        checkbox6.config(state="normal")
        checkbox7.config(state="normal")
        checkbox8.config(state="normal")
        checkbox9.config(state="normal")
        checkbox10.config(state="normal")
        checkbox11.config(state="normal")
        checkbox12.config(state="normal")
        radiobutton1.config(state="normal")
        radiobutton2.config(state="normal")
        boltz_data_frame.config(fg="black")
        ms_word_frame.config(fg="black")
        other_files_frame.config(fg="black")


def energy_temperature_window():
    """Opens a window where the user can modify settings for energy type used in analysis,
    and temperature for Boltzmann-weighted data averaging."""

    def get_parameters():
        """Retrieves newly entered settings and saves these to the settings text file."""
        settings["Temperature (K)"] = entry1.get()
        settings["Boltz energy type"] = var1.get()
        save_new_settings()

    window = Toplevel(root)
    window.iconbitmap(icon_path)
    window.title("Energy & Temperature")
    frame = LabelFrame(window, text="Energy Type Considered")
    frame.pack()
    var1 = StringVar()
    var1.set(settings["Boltz energy type"])
    radiobutton1 = Radiobutton(
        frame,
        text="Electronic energy (E)",
        variable=var1,
        value="Electronic energy",
        anchor="w",
        command=get_parameters,
    )
    radiobutton1.grid(row=2, column=1, sticky="w")
    radiobutton2 = Radiobutton(
        frame,
        text="Gibbs free energy (G)",
        variable=var1,
        value="Gibbs free energy",
        anchor="w",
        command=get_parameters,
    )
    radiobutton2.grid(row=1, column=1, sticky="w")
    space1 = Label(window, text="")
    space1.pack()
    label1 = Label(window, text="Temperature (K) for Boltzmann-weighting")
    label1.pack()
    entry1 = Entry(window, justify="center")
    entry1.insert(END, settings["Temperature (K)"])
    entry1.pack()
    space2 = Label(window, text="")
    space2.pack()
    save_button = Button(window, text=" Save ", command=get_parameters)
    save_button.pack()
    space3 = Label(window, text="")
    space3.pack()


def scaling_factors_settings_window():
    """Opens a window where the user can modify settings for data scaling factors."""

    def get_parameters():
        """Retrieves newly entered settings and saves these to the settings text file."""
        settings["H slope"] = entry1.get()
        settings["H intercept"] = entry2.get()
        settings["C slope"] = entry3.get()
        settings["C intercept"] = entry4.get()
        settings["IR freq scaling factor"] = entry5.get()
        save_new_settings()

    window = Toplevel(root)
    window.iconbitmap(icon_path)
    window.title("Scaling Factors")
    nmr_frame = LabelFrame(window, text="NMR")
    nmr_frame.grid(row=0, column=0)
    nmr_frame.configure(font=("Segoe UI", 9, "bold"))
    ir_frame = LabelFrame(window, text="IR")
    ir_frame.grid(row=1, column=0, pady=15)
    ir_frame.configure(font=("Segoe UI", 9, "bold"))
    label1 = Label(nmr_frame, text="H slope")
    label1.grid(row=0, column=1)
    entry1 = Entry(nmr_frame, justify="center")
    entry1.insert(END, settings["H slope"])
    entry1.grid(row=1, column=1)
    label2 = Label(nmr_frame, text="H intercept")
    label2.grid(row=0, column=3)
    entry2 = Entry(nmr_frame, justify="center")
    entry2.insert(END, settings["H intercept"])
    entry2.grid(row=1, column=3)
    label3 = Label(nmr_frame, text="C slope")
    label3.grid(row=3, column=1)
    entry3 = Entry(nmr_frame, justify="center")
    entry3.insert(END, settings["C slope"])
    entry3.grid(row=4, column=1)
    label4 = Label(nmr_frame, text="C intercept")
    label4.grid(row=3, column=3)
    entry4 = Entry(nmr_frame, justify="center")
    entry4.insert(END, settings["C intercept"])
    entry4.grid(row=4, column=3)
    label5 = Label(ir_frame, text="Frequency slope")
    label5.grid(row=3, column=3)
    entry5 = Entry(ir_frame, justify="center")
    entry5.insert(END, settings["IR freq scaling factor"])
    entry5.grid(row=4, column=3, padx=62)
    save_button = Button(window, text=" Save ", command=get_parameters)
    save_button.grid(row=2, column=0)
    space2 = Label(window, text="")
    space2.grid(row=3, column=0)


def duplicate_conformer_thresholds_window():
    """Opens a window where the user can modify settings for identification of redundant conformers."""

    def get_parameters():
        """Retrieves newly entered settings and saves these to the settings text file."""
        settings["Energy cutoff (kcal/mol)"] = entry1.get()
        settings["MAD cutoff (A)"] = entry2.get()
        settings["Duplicate conformer details"] = var1.get()
        settings["Check for dup confs in XYZ/SDF files"] = var2.get()
        save_new_settings()

    window = Toplevel(root)
    window.title("Redundant conformers")
    window.iconbitmap(icon_path)
    label1 = Label(window, text="Identify redundant conformers by:")
    label1.grid(row=1, column=1)
    label2 = Label(window, text="Energy difference\n(kcal/mol)")
    label2.grid(row=3, column=1)
    entry1 = Entry(window, justify="center")
    entry1.insert(END, settings["Energy cutoff (kcal/mol)"])
    entry1.grid(row=4, column=1)
    label3 = Label(window, text="AND")
    label3.grid(row=5, column=1)
    label4 = Label(window, text="Maximum atom deviation\nbetween aligned geometries\n(MAD, angstroms)")
    label4.grid(row=6, column=1)
    entry2 = Entry(window, justify="center")
    entry2.insert(END, settings["MAD cutoff (A)"])
    entry2.grid(row=7, column=1)
    space2 = Label(window, text="")
    space2.grid(row=8, column=1)
    var1 = BooleanVar()
    var1.set(settings["Duplicate conformer details"])
    checkbox1 = Checkbutton(
        window,
        text="Save energy differences and MADs of\nredundant conformers to a .txt file",
        variable=var1,
        anchor="w",
        command=get_parameters,
    )
    checkbox1.grid(row=11, column=1, sticky="w")
    var2 = BooleanVar()
    var2.set(settings["Check for dup confs in XYZ/SDF files"])
    checkbox2 = Checkbutton(
        window,
        text="Look for redundant conformers from\n.xyz/.sdf files too (MAD only)",
        variable=var2,
        anchor="w",
        command=get_parameters,
    )
    checkbox2.grid(row=12, column=1, sticky="w")
    space3 = Label(window, text="")
    space3.grid(row=10, column=1)
    save_button = Button(window, text=" Save ", command=get_parameters)
    save_button.grid(row=9, column=1)


def input_file_writer_window(
    settings,
    data,
    initial_filename,
    multi_out_file_flag,
    conformer_renumber_flag,
    keep_imag_freq_flag,
    out_file_content_types,
):
    """Opens a window with an interactive template for creating quantum chemistry calculation input files."""

    global input_file_texts
    global inp_window_user_decision
    inp_window_user_decision = ""
    input_file_texts = settings["Input File Texts"]
    calc_sections = []
    all_entries = []
    global add_calc_count
    add_calc_count = 0
    global message_height
    message_height = 0
    global compound_name
    var1 = BooleanVar()
    var1.set(settings["Input File Conformers Together"])
    inp_writer_window = Toplevel(root)
    inp_writer_window.title("Input File Template")
    inp_writer_window.iconbitmap(icon_path)
    inp_writer_window.maxsize(664, 1000)
    # Set this below setting to false, to start with.

    settings["Skip excluding duplicate conformers from input files made from output files"] = False

    def get_parameters():
        """Retrieves newly entered settings and saves these to the settings text file."""
        settings["Input File Conformers Together"] = var1.get()
        global input_file_texts
        settings["Input File Texts"] = input_file_texts
        global compound_name
        compound_name = filename_entry.get()
        if len(data) != 6:  # This means data is from .out/.log files
            settings["Geometry reopt relative energy threshold"] = entry0_var.get()
            settings["Geometry reopt relative energy unit"] = var_unit.get()
            settings["Geometry reopt"] = var_reopt.get()
            settings["Skip excluding duplicate conformers from input files made from output files"] = (
                var_not_excl_dup_confs.get()
            )
            if var_reopt.get() is False:
                entry0.config(state="disabled")
                radiobutton1.config(state="disabled")
                radiobutton2.config(state="disabled")
                reopt_tick.config(fg="grey")
            else:
                entry0.config(state="normal")
                radiobutton1.config(state="normal")
                radiobutton2.config(state="normal")
                reopt_tick.config(fg="black")
        update_widgets()
        save_new_settings()

    def create_input_file(all_entries):
        """Retrieves template data, closes window and tells main() that user has selected to create input file(s)
        with this data."""

        global inp_window_user_decision
        inp_window_user_decision = "Create input file."
        global input_file_texts
        input_file_texts = []
        texts = []
        texts.append(spacer_entry.get("1.0", "end-1c"))
        for object in all_entries:
            first_part = ""
            last_part = ""
            if str(object).endswith("!labelframe.!text"):
                first_part = object.get("1.0", "end-1c")
            elif str(object).endswith("!labelframe.!text2"):
                last_part = object.get("1.0", "end-1c")
                texts.append(first_part + last_part)
            else:
                texts.append(object.get("1.0", "end-1c"))
        input_file_texts = texts
        get_parameters()
        inp_writer_window.destroy()
        return

    def update_widgets():
        """Updates widgets in input file template window with currently entered values."""

        global message_height
        try:
            filename_label.configure(text="Compound name: ")  # This throws an error on start-up but is fine afterwards
        except:
            return
        if var1.get() is True:
            filename_label.configure(text="Compound name: ")
            if len(data) != 6 and var_reopt.get() is True:  # This means data is from output files and an energy window
                # is selected

                filename_label2_text = filename_entry.get() + "_" + str(entry0.get()) + str(var_unit.get()) + ".inp"
            else:
                filename_label2_text = filename_entry.get() + ".inp"
            filename_label2.configure(text=filename_label2_text)
            filename_label_text = "Filename: "
            filename_label3.configure(text=filename_label_text)
            Button2.configure(text=" Save template & create input file ")
        elif var1.get() is False:
            filename_label_text = "Filenames: "
            if len(data) != 6 and var_reopt.get() is True:  # This means data is from output files and an energy window
                # is selected

                filename_label2_text = (
                    filename_entry.get()
                    + "_"
                    + str(entry0.get())
                    + str(var_unit.get())
                    + "_conf-⫷⫷⫷CONFORMER NUMBER⫸⫸⫸.inp"
                )
            else:
                filename_label2_text = filename_entry.get() + "_conf-⫷⫷⫷CONFORMER NUMBER⫸⫸⫸.inp"
            filename_label2.configure(text=filename_label2_text)
            filename_label3.configure(text=filename_label_text)
            Button2.configure(text=" Save template & create input files ")
        if len(data) != 6:  # This means data is from .out/.log files
            # Find number of confs below energy threshold

            rel_energy_threshold = float(entry0_var.get())
            if (
                settings["Geometry reopt relative energy unit"] == "kJmol-1"
                and settings["Relative energy unit"] == "kcal/mol"
            ):
                rel_energy_threshold = rel_energy_threshold / 4.184
            elif (
                settings["Geometry reopt relative energy unit"] == "kcalmol-1"
                and settings["Relative energy unit"] == "kJ/mol"
            ):
                rel_energy_threshold = rel_energy_threshold * 4.184
            number_confs_above_threshold = str(sum(i > rel_energy_threshold for i in data[1]))
            if number_confs_above_threshold == "1":
                reopt_tick.configure(
                    text=number_confs_above_threshold + " unique conformer with relative energy above this\nvalue will "
                    "be 𝐞𝐱𝐜𝐥𝐮𝐝𝐞𝐝 from the calculations below."
                )
            else:
                reopt_tick.configure(
                    text=number_confs_above_threshold + " unique conformers with relative energies above this\nvalue "
                    "will be 𝐞𝐱𝐜𝐥𝐮𝐝𝐞𝐝 from the calculations below."
                )
            if "program_error_result" in globals():
                if program_error_result == "proceed without exclusion":
                    reopt_tick.configure(
                        text=number_confs_above_threshold + " conformers (excluding conformers with "
                        "imaginary\nfrequencies) with relative energies above this "
                        "value\nwill be excluded from the calculations below."
                    )
            reopt_tick.update()
            # Update warning to user if subsequent input files would be backwards-incompatible, and propose solutions

            incompatibility_frame.grid_forget()
            no_dup_conf_exclusion_box.pack_forget()
            message_height = 0
            if keep_imag_freq_flag != "proceed without exclusion":
                if not data[50] or var1.get() is True:
                    if multi_out_file_flag == "Multiple output files" and var1.get() is False:
                        incompatibility_frame.grid(row=199, column=1)
                        incompatibility_label.pack()
                        incompatibility_label.configure(font=("Segoe UI", 9))
                        incompatibility_text = (
                            "CAUTION: Backwards-incompatible!\nYour current .out\.log filenames do not have"
                            " conformer suffixes (e.g. conf-XX). This means that calculated "
                            "results\nfrom the new input files will be incompatible for combined "
                            "analysis with your current .out\.log file(s).\nThis incompatibility may be "
                            "fine (e.g. if doing another geometry optimization). "
                        )
                        incompatibility_label.configure(text=incompatibility_text)
                        message_height = 95
                    elif multi_out_file_flag == "" and var1.get() is False:
                        incompatibility_frame.grid(row=199, column=1)
                        incompatibility_label.pack()
                        incompatibility_label.configure(font=("Segoe UI", 9))
                        incompatibility_text = (
                            "CAUTION: Backwards-incompatible!\nYou have selected to have conformers in"
                            " separate input files, but your current .out\.log file has all conformers in\n"
                            "the same file. Calculated results from these new input files will be "
                            "incompatible for combined analysis with the\ncurrent .out\.log file(s). This "
                            "is because conformer data can not be matched across these new and old "
                            "files with\ndifferent conformer arrangements during data analysis. This "
                            "incompatibility may be fine (e.g. if doing another\ngeometry "
                            "optimization), or can be resolved now by selecting to have all conformers"
                            " in the same input file. "
                        )
                        incompatibility_label.configure(text=incompatibility_text)
                        message_height = 95
                    elif multi_out_file_flag and var1.get() is True and len(set(out_file_content_types)) == 1:
                        incompatibility_frame.grid(row=199, column=1)
                        incompatibility_label.pack()
                        incompatibility_label.configure(font=("Segoe UI", 9))
                        incompatibility_text = (
                            "CAUTION: Backwards-incompatible!\nYou have selected to have all "
                            "conformers in the same input file, but your current .out\.log files have "
                            "conformers\nin separate files. Calculated results from this new input "
                            "file will be incompatible for combined analysis later with\nthe current "
                            ".out\.log files. This is because conformer data can not be matched across "
                            "these new and old files with\ndifferent conformer arrangements during "
                            "data analysis. This incompatibility may be fine (e.g. if doing "
                            "another\ngeometry optimization), or can be resolved now by deselecting to"
                            " have all conformers in the same input file. "
                        )
                        incompatibility_label.configure(text=incompatibility_text)
                        message_height = 95
                    elif (
                        data[32] or settings["Geometry reopt"] is True and int(number_confs_above_threshold) > 0
                    ):  # If dup or high-energy conformers are going to be
                        # removed (i.e., conf ordering changed)

                        incompatibility_frame.grid(row=199, column=1)
                        incompatibility_label.pack()
                        incompatibility_text = "CAUTION: Backwards-incompatible!\nExcluding "
                        a = ""
                        if data[32]:  # Dup confs to be removed
                            if int(data[32][0]) == 1:
                                a += str(data[32][0]) + " redundant conformer "
                            elif int(data[32][0]) > 1:
                                a += str(data[32][0]) + " redundant conformers "
                        if (
                            data[32] and settings["Geometry reopt"] is True and int(number_confs_above_threshold) > 0
                        ):  # both
                            a += "and "
                        if (
                            settings["Geometry reopt"] is True and int(number_confs_above_threshold) > 0
                        ):  # High-energy confs to be removed
                            if int(number_confs_above_threshold) == 1:
                                a += number_confs_above_threshold + " high-energy conformer "
                                incompatibility_text += (
                                    a + "from this new input file will make the resulting\n.out\.log file "
                                    "incompatible for combined data analysis with the current "
                                    ".out\.log file(s). This is because excluding this\n "
                                    "conformer changes the ordering of conformers in the input "
                                    "file, preventing the correct matching of conformer data\n"
                                    "across output files during data analysis in SpectroIBIS. This "
                                    "incompatibility may be fine (e.g. if doing another geometry\n"
                                    "optimization), or can be resolved now by deselecting "
                                    "the conformer energy window above and ticking the box below. "
                                )
                            elif int(number_confs_above_threshold) > 1:
                                a += number_confs_above_threshold + " high-energy conformers "
                                incompatibility_text += (
                                    a + "from this new input file will make the resulting\n.out\.log file "
                                    "incompatible for combined data analysis with the current "
                                    ".out\.log file(s). This is because excluding these\n "
                                    "conformers changes the ordering of conformers in the input "
                                    "file, preventing the correct matching of conformer data\n"
                                    "across output files during data analysis in SpectroIBIS. This "
                                    "incompatibility may be fine (e.g. if doing another geometry\n"
                                    "optimization), or can be resolved now by deselecting "
                                    "the conformer energy window above and ticking the box below. "
                                )
                        else:
                            incompatibility_text += (
                                a + "from this new input file will make the resulting .out\.log\nfile "
                                "incompatible for combined data analysis with the current .out\.log "
                                "file(s). This is because excluding\nthese conformer(s) changes "
                                "the ordering of conformers in the input file, preventing the "
                                "correct matching of\nconformer data across output files during data "
                                "analysis in SpectroIBIS. This incompatibility may be fine (e.g."
                                "\nif doing another geometry optimization), or can be resolved now"
                                " by ticking the box below. "
                            )
                        incompatibility_label.configure(text=incompatibility_text)
                        no_dup_conf_exclusion_box.pack()
                        message_height = 145
            # Tell user if new conformer numbers will be generated.

            if conformer_renumber_flag is True:
                conformer_renumber_frame.grid(row=198, column=1)
                conformer_renumber_label.pack()
                conformer_renumber_label.configure(font=("Segoe UI", 9))
                conformer_renumber_text = (
                    "\nCAUTION: Conformer numbers may differ from current .out\.log file(s).\nBe careful if using checkpoint "
                    "files (Gaussian) or .xyz file geometry inputs (ORCA). This happens because no\nconformer numbers "
                    "could be extracted from these .out\.log files, so new conformer numbers must be arbitrarily assigned. "
                )
                conformer_renumber_label.configure(text=conformer_renumber_text)
                message_height += 75
            # Update window height

            inp_window_height = 440 + 192 * (int(add_calc_count / 5) - 1)
            inp_window_height += 96 + message_height
            inp_window_geom = "703x" + str(inp_window_height)
            inp_writer_window.geometry(inp_window_geom)

    def remove_calc(calc_sections, all_entries):
        """Removes one calculation from input file template."""

        global add_calc_count
        global message_height
        if add_calc_count <= 5:
            return
        else:
            calc_sections[-1][0].grid_forget()
            del calc_sections[-1]
            del all_entries[-1]
            add_calc_count -= 5
        inp_window_height = 440 + 192 * (int(add_calc_count / 5) - 1)
        if len(data) != 6:  # This means data is from .out/.log files
            inp_window_height += 96 + message_height
        inp_window_geom = "703x" + str(inp_window_height)
        inp_writer_window.geometry(inp_window_geom)

    def add_calc():
        """Adds one calculation to input file template."""

        global add_calc_count
        global message_height
        add_calc_count += 5
        frame = LabelFrame(inp_inner_frame, text="Calculation #" + str(int(add_calc_count / 5)), padx=10, pady=4)
        frame.configure(font=("Segoe UI", 9, "bold"))
        frame.grid(row=(11 + add_calc_count), column=1, padx=10, pady=10)
        if add_calc_count == 5:
            entry_height = 9
            main_text = settings["Input File Texts"][int(add_calc_count / 5)]
        elif add_calc_count > 5:
            entry_height = 8
            spacer_label = Label(frame, text="Multiple job text will be pasted here.", justify="left")
            spacer_label.grid(row=(12 + add_calc_count), column=0, sticky="w")
            if int(add_calc_count / 5) <= (len(settings["Input File Texts"]) - 1):
                main_text = settings["Input File Texts"][int(add_calc_count / 5)]
            else:
                main_text = ""
        entry_main = Text(frame, width=100, height=entry_height, undo=True)
        entry_main.insert(END, main_text)
        entry_main.grid(row=(14 + add_calc_count), column=0)
        entry_main.configure(font=("Segoe UI", 9))
        all_entries.append(entry_main)
        calc_section = [frame]
        calc_sections.append(calc_section)
        inp_window_height = 440 + 192 * (int(add_calc_count / 5) - 1)
        if len(data) != 6:  # This means data is from .out/.log files
            inp_window_height += 96 + message_height
        inp_window_geom = "703x" + str(inp_window_height)
        inp_writer_window.geometry(inp_window_geom)

    inp_outer_frame = Frame(inp_writer_window)
    inp_outer_frame.pack(fill="both", expand=1)
    inp_canvas = Canvas(inp_outer_frame)
    inp_canvas.pack(side="left", fill="both", expand=1)
    inp_scrollbar = Scrollbar(inp_outer_frame, orient="vertical", command=inp_canvas.yview)
    inp_scrollbar.pack(side="right", fill="y")
    inp_canvas.configure(yscrollcommand=inp_scrollbar.set)
    inp_canvas.bind("<Configure>", lambda e: inp_canvas.configure(scrollregion=inp_canvas.bbox("all")))
    inp_inner_frame = Frame(inp_canvas)
    inp_canvas.create_window((0, 0), window=inp_inner_frame, anchor="nw")
    Button2 = Button(inp_inner_frame, text=" Create input file(s) ", command=lambda: create_input_file(all_entries))

    def callback(entry_var):
        """Triggers widgets to update with newly enetered value(s)."""

        update_widgets()

    def test_val(entered_text, action_type):
        """Validates inputted character is a digit or dot point."""

        if action_type == "1":  # insertion of a character
            if entered_text.endswith("."):
                return True
            if not entered_text[-1].isdigit():
                return False
        return True

    # Add option to exclude high-energy conformers

    if len(data) != 6:  # This means data is from .out/.log files
        intro_message = str(len(data[0])) + " total"
        global program_error_result
        if "program_error_result" in globals():
            if program_error_result == "proceed without exclusion":  # This means imag freq confs were detected and user
                # has selected to not exclude them.

                number_imag_freq_confs = data[53]
                intro_message = str(len(data[0]) + number_imag_freq_confs) + " total"
        if data[32]:  # dup confs present
            if data[32][0] == "1":
                intro_message += ", excluding " + str(data[32][0]) + " redundant conformer"
            else:
                intro_message += ", excluding " + str(data[32][0]) + " redundant conformers"
        if keep_imag_freq_flag == "proceed without exclusion":
            if data[32]:  # dup confs present
                intro_message = str(len(data[0]) + number_imag_freq_confs + int(data[32][0])) + " total"
            else:
                intro_message = str(len(data[0]) + number_imag_freq_confs) + " total"
        var_reopt = BooleanVar()
        var_reopt.set(settings["Geometry reopt"])
        var_unit = StringVar()
        var_unit.set(settings["Geometry reopt relative energy unit"])
        var_not_excl_dup_confs = BooleanVar()
        reopt_frame = LabelFrame(inp_inner_frame, text="Conformer Energy Window")
        reopt_frame.configure(font=("Segoe UI", 9, "bold"))
        reopt_frame.grid(row=1, column=1, pady=5)
        reopt_frame2 = LabelFrame(reopt_frame, borderwidth=0, highlightthickness=0)
        reopt_frame2.pack()
        entry0_var = StringVar()
        entry0_var.trace("w", lambda name, index, mode, sv=entry0_var: callback(entry0_var))
        entry0 = Entry(reopt_frame2, validate="key", textvariable=entry0_var, justify="center", width=10)
        entry0["validatecommand"] = (entry0.register(test_val), "%P", "%d")
        entry0.insert(END, settings["Geometry reopt relative energy threshold"])
        entry0.grid(row=0, column=1)
        radiobutton1 = Radiobutton(
            reopt_frame2, text="kJ/mol", variable=var_unit, value="kJmol-1", anchor="w", command=get_parameters
        )
        radiobutton1.grid(row=0, column=2, sticky="w")
        radiobutton2 = Radiobutton(
            reopt_frame2, text="kcal/mol", variable=var_unit, value="kcalmol-1", anchor="w", command=get_parameters
        )
        radiobutton2.grid(row=0, column=3, sticky="w")
        reopt_frame3 = LabelFrame(reopt_frame, borderwidth=0, highlightthickness=0)
        reopt_frame3.pack()
        # Find number of confs below energy threshold

        rel_energy_threshold = float(settings["Geometry reopt relative energy threshold"])
        if settings["Geometry reopt relative energy unit"] == "kJmol-1" and settings[
            "Relative energy unit"] == "kcal/mol":
            rel_energy_threshold = rel_energy_threshold / 4.184
        elif (
                settings["Geometry reopt relative energy unit"] == "kcalmol-1" and settings[
            "Relative energy unit"] == "kJ/mol"
        ):
            rel_energy_threshold = rel_energy_threshold * 4.184
        number_confs_above_threshold = str(sum(i > rel_energy_threshold for i in data[1]))
        reopt_tick = Checkbutton(
            reopt_frame3,
            variable=var_reopt,
            anchor="w",
            text=number_confs_above_threshold + " conformers with relative energies above this "
                                                "value\nwill be excluded from the calculations "
                                                "below.",
            command=get_parameters,
        )
        reopt_tick.grid(row=0, column=0)
        # If user isn't excluding non-converged conformers, they wouldn't want to exclude duplicates or
        # high-energy conformers either.

        if keep_imag_freq_flag == "proceed without exclusion":
            var_reopt.set(False)
            var_not_excl_dup_confs.set(True)
            reopt_frame.grid_forget()
        # Update relevant widgets

        if var_reopt.get() is False:
            entry0.config(state="disabled")
            radiobutton1.config(state="disabled")
            radiobutton2.config(state="disabled")
            reopt_tick.config(fg="grey")
        else:
            entry0.config(state="normal")
            radiobutton1.config(state="normal")
            radiobutton2.config(state="normal")
            reopt_tick.config(fg="black")
        # Define a frame to be packed later, to warn user if subsequent input files would be backwards-incompatible
        # and propose solutions

        incompatibility_frame = LabelFrame(inp_inner_frame, borderwidth=0, highlightthickness=0)
        incompatibility_text = ""
        incompatibility_label = Label(incompatibility_frame, text=incompatibility_text, fg="blue")
        no_dup_conf_exclusion_box = Checkbutton(
            incompatibility_frame,
            text="Don't exclude redundant or high-energy conformers from this "
                 "input\nfile to maintain analysis-compatibility with current "
                 ".out/.log file(s).\n(Redundant conformers will be excluded during "
                 "data analysis later.)",
            variable=var_not_excl_dup_confs,
            command=get_parameters,
        )
        # Define a frame to be packed later, to warn user if conformers will be renumbered.

        conformer_renumber_frame = LabelFrame(inp_inner_frame, borderwidth=0, highlightthickness=0)
        conformer_renumber_text = ""
        conformer_renumber_label = Label(conformer_renumber_frame, text=conformer_renumber_text)
    if len(data) == 6:  # This means data is from .xyz/.sdf files
        intro_message = str(len(data[0][0])) + " total"
        if data[5]:  # dup confs present
            if data[5][0] == "1":
                intro_message += ", excluding " + str(data[5][0]) + " redundant conformer"
            else:
                intro_message += ", excluding " + str(data[5][0]) + " redundant conformers"
    intro_label = Label(
        inp_inner_frame,
        text="Please enter input file details into the template below.\nThis template will be used "
             "for each unique conformer (" + intro_message + ").\nCompound name, relevant conformer "
                                                             "numbers and atom coordinates will be "
                                                             'pasted wherever the terms\n"⫷⫷⫷COMPOUND '
                                                             'NAME⫸⫸⫸", "⫷⫷⫷CONFORMER NUMBER⫸⫸⫸"and '
                                                             '"⫷⫷⫷ATOM COORDINATES⫸⫸⫸"\nare found, '
                                                             "respectively. These terms can be added "
                                                             "and removed as desired.",
    )
    intro_label.grid(row=0, column=1)
    intro_label.configure(font=("Segoe UI", 9))
    first_entries_frame = LabelFrame(inp_inner_frame, borderwidth=0, highlightthickness=0)
    first_entries_frame.grid(row=3, column=1, pady=5)
    second_entries_frame = LabelFrame(inp_inner_frame, borderwidth=0, highlightthickness=0, pady=5)
    second_entries_frame.grid(row=5, column=1)
    if var1.get() is True:
        filename_label_text = "Filename: "
        filename_label2_text = ".inp"
    else:
        filename_label_text = "Filenames: "
        filename_label2_text = "_conf-⫷⫷⫷CONFORMER NUMBER⫸⫸⫸.inp"
    filename_label = Label(first_entries_frame, text=filename_label_text)
    filename_label.grid(row=0, column=0)
    filename_label.configure(font=("Segoe UI", 9, "bold"))
    name_var = StringVar(value=initial_filename)
    name_var.trace("w", lambda name, index, mode, sv=name_var: callback(name_var))
    filename_entry = Entry(first_entries_frame, textvariable=name_var, width=51)
    filename_entry.grid(row=0, column=1)
    filename_entry.configure(font=("Segoe UI", 9))
    filename_frame = LabelFrame(inp_inner_frame, borderwidth=0, highlightthickness=0, bg="red")
    filename_frame.grid(row=4, column=1)
    filename_label2 = Label(filename_frame, text=filename_label2_text)
    filename_label2.grid(row=0, column=1)
    filename_label3 = Label(filename_frame, text="Filename: ")
    filename_label3.grid(row=0, column=0)
    filename_label3.configure(font=("Segoe UI", 9, "bold"))
    update_widgets()
    spacer_label = Label(second_entries_frame, text="Multiple job text: \n(when applicable)")
    spacer_label.grid(row=1, column=1)
    spacer_label.configure(font=("Segoe UI", 9, "bold"))
    spacer_entry = Text(second_entries_frame, width=20, height=2, undo=True)
    spacer_entry.insert(END, settings["Input File Texts"][0])
    spacer_entry.grid(row=1, column=2)
    spacer_entry.configure(font=("Segoe UI", 9))
    if all_entries:
        if all_entries[0] != spacer_entry:
            all_entries.append(spacer_entry)
    checkbox_conformers_together = Checkbutton(
        second_entries_frame,
        text="All conformers in same input file                         ",
        variable=var1,
        command=update_widgets,
    )
    checkbox_conformers_together.grid(row=1, column=0)
    checkbox_conformers_together.configure(font=("Segoe UI", 9, "bold"))
    # Open a number of calc boxes that matches what was last saved in the settings file.

    for i in range(0, (len(settings["Input File Texts"]) - 1)):
        add_calc()
    button_frame = LabelFrame(inp_inner_frame, borderwidth=0, highlightthickness=0)
    button_frame.grid(row=197, column=1)
    add_calc_button = Button(button_frame, text=" Add another calculation ", command=add_calc)
    add_calc_button.grid(row=0, column=0, padx=5)
    remove_calc_button = Button(
        button_frame, text=" Remove last calculation ", command=lambda: remove_calc(calc_sections, all_entries)
    )
    remove_calc_button.grid(row=0, column=1, padx=5)
    Button2.grid(row=202, column=1, pady=10)
    Button2.configure(font=("Segoe UI", 9, "bold"))
    inp_writer_window.wait_window(inp_writer_window)
    try:
        compound_name
    except NameError:
        compound_name = ""
    return inp_window_user_decision, compound_name


# Create an info menu


info_menu = Menu(menu_bar, tearoff=0)
menu_bar.add_cascade(label="Info", menu=info_menu)
info_menu.add_command(label="User Manual (PDF)", command=open_manual)
info_menu.add_command(label="About", command=open_about_window)

# Create a settings menu


settings_menu = Menu(menu_bar, tearoff=0)
menu_bar.add_cascade(label="Settings", menu=settings_menu)
settings_menu.add_command(label="Output", command=output_window)
settings_menu.add_command(label="Redundant Conformers", command=duplicate_conformer_thresholds_window)
settings_menu.add_command(label="Scaling Factors", command=scaling_factors_settings_window)
settings_menu.add_command(label="Energy & Temperature", command=energy_temperature_window)


# Create a status bar


status_text = ""
if settings["Mode"] == "Create input files":
    status_text = "Note: SpectroIBIS is currently set to create input file(s) from output files."
status_bar = Label(root, text=status_text)
status_bar.pack()

# Create a second status bar


status_text2 = ""
if settings["Mode"] == "Create input files":
    status_text2 = "This can be changed in Settings --> Output."  # Warn user that SpectroIBIS is set to create input
    # files from .out/.log files
status_bar2 = Label(root, text=status_text2)
status_bar2.pack()

# Make window drag-and-droppable


root.drop_target_register(DND_FILES)
root.dnd_bind("<<Drop>>", on_drop)

# Start the tkinter main loop


root.mainloop()
