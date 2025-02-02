import os
import re
from os import listdir
from os.path import isfile, join


def process_string_from_header(input_string):
    # Define the function to replace based on the condition
    def replace_match(match):
        # If it's an integer, remove the dot and integer
        if match.group(2).isdigit():
            return match.group(1)
        # If it's not an integer, replace the dot with a hyphen
        return match.group(1) + "-" + match.group(2)

    # Use regex to find patterns with a dot followed by any characters
    result = re.sub(r'(\w+)\.(\w+)', replace_match, input_string)
    return result

def multiline_fasta_check(file):
    with open(file, "r") as f:
        lines = f.readlines()
    number_of_inputs = sum([1 for line in lines if ">" in line])
    return number_of_inputs != 1


def multiline_fasta_handle(file):
    base_name = str(os.path.basename(file).split(".")[0])
    try:
        os.mkdir(base_name)
    except OSError:
        pass

    cmd = f"cat {file}"
    cmd += " | awk '{ if (substr($0, 1, 1)==\">\") {"
    cmd += "filename=(\"{}/\"".format(base_name)
    cmd += "substr($0,2)\".fa\")} "
    cmd += f"print $0 > filename "
    cmd += "}'"

    os.system(cmd)

    return base_name


def multiline_fasta_handle_python(file, flag_ncbi_formatting=False):
    base_name = str(os.path.basename(file).split(".")[0])
    try:
        os.mkdir(base_name)
    except OSError:
        pass

    with open(file, "r") as f:
        lines = f.readlines()

    headers = []
    dna_sequences = []

    dna_sequence = ''
    for line in lines:
        if line:
            if ">" in line:
                if dna_sequence:
                    dna_sequences.append(dna_sequence)
                    dna_sequence = ''
                headers.append(line)
            else:
                dna_sequence += line.strip()

    if dna_sequence:
        dna_sequences.append(dna_sequence)

    if flag_ncbi_formatting:
        for header, dna_sequence in zip(headers, dna_sequences):
            new_header = header.split(" ")[0]
            new_header = process_string_from_header(new_header)
            file_name = new_header.split(">")[1].replace(",", "-") \
                                  .replace(".", "-").replace(" ", "_").replace("|", "-") + ".fa"
            with open(os.path.join(base_name, file_name), "w") as f:
                f.writelines(new_header)
                f.write("\n")
                f.writelines(dna_sequence)
    else:
        for header, dna_sequence in zip(headers, dna_sequences):
            file_name = header.strip().split(">")[1].replace(",", "_")\
                            .replace(".", "_").replace(" ", "_").replace("|", "_") + ".fa"
            with open(os.path.join(base_name, file_name), "w") as f:
                f.write(header)
                f.write(dna_sequence)

    return base_name


def folder_of_multifasta_handle(folder_multifasta):
    list_files = [f for f in listdir(folder_multifasta) if isfile(join(folder_multifasta, f))]
    all_lines_in_files = []
    for file in list_files:
        with open(os.path.join(folder_multifasta, file), "r") as f:
            lines = f.readlines()
        all_lines_in_files.append(lines)
    with open("multifasta_folder.fa", "w") as f:
        for lines in all_lines_in_files:
            for line in lines:
                f.write(line)

    multiline_fasta_handle_python("multifasta_folder.fa")
    return "multifasta_folder"