import argparse
import sys
from os import listdir
from os.path import isfile, join
import os
from time import time
import multiprocessing
import subprocess

sys.path.insert(0, 'components/')
from predictor import Predictor
from machine_learning import ClassifierWrapper
import shutil
from time import time
import math


parser = argparse.ArgumentParser(description='Run Identifier')
parser.add_argument('--input_folder', nargs='*', type=str, default=None,
                    help='input folder (default: None)')

parser.add_argument('--file', type=str, default=None,
                    help='input file (default: None)')

parser.add_argument('--model', type=str, default="8",
                    help='model_to_use (default: 8)')

parser.add_argument('--result_folder', type=str, default="Results",
                    help='folder with the result (default: Results)')

parser.add_argument('--pickle_report', type=str, default=None,
                    help='pickled report file (default: None)')

parser.add_argument('--cas', type=bool, default=False,
                    help='cas genes computation (default: False)')

parser.add_argument('--is_element', type=bool, default=False,
                    help='is element computation (default: False)')

parser.add_argument('--parallel', type=str, default=True,
                    help='parallel computations (default: True)')

parser.add_argument('--cpu', type=str, default="ALL",
                    help='parallel computations (default: ALL)')

parser.add_argument('--fast_run', type=str, default=False,
                    help='fast run option (default: False)')

parser.add_argument('--degenerated', type=bool, default=True,
                    help='degenerated_repeat_computation (default: False)')

parser.add_argument('--min_len_rep', type=int, default=21,
                    help='min avg. length of the repeats (default: 21)')

parser.add_argument('--max_len_rep', type=int, default=55,
                    help='max avg. length of the repeats (default: 55)')

parser.add_argument('--min_len_spacer', type=int, default=18,
                    help='min avg. length of spacers (default: 18)')

parser.add_argument('--max_len_spacer', type=int, default=78,
                    help='max avg. length of spacers (default: 78)')

parser.add_argument('--min_repeats', type=int, default=3,
                    help='min number of repeats (default: 3)')

parser.add_argument('--log_file', type=str, default=None,
                    help='log file (default: None)')

parser.add_argument('--enhancement_max_min', type=bool, default=True,
                    help='enhancement with filter (default: True)')

parser.add_argument('--enhancement_start_end', type=bool, default=True,
                    help='enhancement with start end omitting (default: True)')


args = parser.parse_args()
complete_path_folder_info = args.input_folder
complete_path_folder = complete_path_folder_info[0]

if len(complete_path_folder_info) == 3:
    chunk_number = int(complete_path_folder_info[1])
    number_of_chunks = int(complete_path_folder_info[2])
else:
    chunk_number = None
    number_of_chunks = None

complete_path_file = args.file
folder_result = args.result_folder
report_pickle = args.pickle_report
list_models = ["8", "9", "10"] if args.model == "ALL" else [args.model]


#flag_parallel = args.parallel
#flag_cpu = args.cpu
#cas_flag = args.cas
#is_flag = args.is_element
#degenerated_flag = args.degenerated
#flag_fast_run = args.fast_run
flag_enhancement_max_min = args.enhancement_max_min
flag_enhancement_start_end = args.enhancement_start_end

flag_parallel = False if (args.parallel in ["False", False]) else True
flag_cpu = args.cpu
flag_fast_run = False if (args.fast_run in ["False", False]) else True

cas_flag = False if (args.cas in ["False", False]) else True
is_flag = False if (args.is_element in ["False", False]) else True
degenerated_flag = False if (args.degenerated in ["False", False]) else True

min_rep = args.min_len_rep
max_rep = args.max_len_rep
max_spacer = args.max_len_spacer
min_spacer = args.min_len_spacer
min_repeats = args.min_repeats

parameters = {
    "param_min_avg_repeat_length": min_rep,
    "param_max_avg_repeat_length": max_rep,
    "param_max_avg_spacer_length": max_spacer,
    "param_min_avg_spacer_length": min_spacer,
    "param_min_repeats": min_repeats
}

log_file = args.log_file


ALL_FEATURES = ['repeat_len', 'number_repeats', 'repeat_similarity',
                'at_richness', 'avg_spacer_len', 'spacer_similarity',
                'number_mismatches', 'spacer_evenness', 'mfe_score',
                'orf_score', 'hmmr_score', 'blast_score_1', 'blast_score_2',
                'eden_score']

best_combinations = {
    "8": (2, 4, 5, 6, 7, 8, 9, 11),
    "9": (1, 2, 4, 5, 7, 8, 9, 10, 12),
    "10": (0, 2, 3, 4, 5, 6, 7, 10, 11, 12)
}

feature_list = ['.'.join([ALL_FEATURES[i] for i in best_combinations[model]]) for model in list_models]
list_ml_classifiers = [ClassifierWrapper(classifier_type=None,
                                         load_option="trained_models/extra_trees/extra_trees_subset{}features.pkl".
                                                     format(model))
                       for model in list_models]


def run_over_folder_of_files(folder, result_folder, chunk_number, number_of_chunks, report_pickle_folder=report_pickle):
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    files = sorted(files)

    if number_of_chunks:
        chunk_size = math.ceil(len(files) / number_of_chunks)
        chunk_start = (chunk_number - 1) * chunk_size
        chunk_end = chunk_number * chunk_size
        chunk = files[chunk_start:chunk_end]
        print(chunk_start)
        print(chunk_end)
    else:
        chunk = files

    for index, file in enumerate(chunk, 1):
        print("\n\n\n\t\t\t\tExecuting file {} out of {} ({})\n\n\n".format(index, len(chunk), file))
        pr = Predictor(result_folder_path="{}/".format(result_folder),
                       file_path="{}/{}".format(folder, file),
                       eden_classifier=None, list_ml_classifiers=list_ml_classifiers,
                       list_features=feature_list, flag_is=is_flag, flag_cas=cas_flag,
                       flag_degenerated=degenerated_flag,
                       flag_parallel=flag_parallel,
                       flag_cpu=flag_cpu,
                       flag_fast_run=flag_fast_run,
                       flag_enhancement_max_min=flag_enhancement_max_min,
                       flag_enhancement_start_end=flag_enhancement_start_end,
                       parameters=parameters,
                       log_file=log_file
                       )

        if report_pickle_folder:
            pr.report_into_separate_file(report_pickle_folder)


def run_over_one_file(file, result_folder, report_pickle_folder=report_pickle):
    print("\n\n\n\t\t\t\tExecuting file {}\n\n\n".format(file))
    pr = Predictor(result_folder_path="{}/".format(result_folder),
                   file_path="{}".format(file),
                   eden_classifier=None, list_ml_classifiers=list_ml_classifiers,
                   list_features=feature_list, flag_is=is_flag, flag_cas=cas_flag,
                   flag_degenerated=degenerated_flag,
                   flag_parallel=flag_parallel,
                   flag_cpu=flag_cpu,
                   flag_fast_run=flag_fast_run,
                   parameters=parameters,
                   flag_enhancement_max_min=flag_enhancement_max_min,
                   flag_enhancement_start_end=flag_enhancement_start_end,
                   log_file=log_file)

    if report_pickle_folder:
        pr.report_into_separate_file(report_pickle_folder)


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
    cmd += "substr($0,2) \".fa\")} "
    cmd += f"print $0 > filename "
    cmd += "}'"

    return base_name


if __name__ == "__main__":
    start_time = time()
    if complete_path_file:
        if multiline_fasta_check(complete_path_file):
            folder_multifasta = multiline_fasta_handle(complete_path_file)
            run_over_folder_of_files(folder_multifasta, folder_result,
                                     chunk_number=chunk_number, number_of_chunks=number_of_chunks)
            shutil.rmtree(folder_multifasta)
        else:
            run_over_one_file(complete_path_file, folder_result)
    else:
        run_over_folder_of_files(complete_path_folder, folder_result,
                                 chunk_number=chunk_number, number_of_chunks=number_of_chunks)

    end_time = time()
    print("Elapsed time: ", end_time-start_time)
