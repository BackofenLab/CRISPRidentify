import argparse
import math
import shutil
import warnings
import os

from pathlib import Path
from os import listdir
from os.path import isfile, join
from time import time

from components.pipeline import Pipeline
from components.components_ml import ClassifierWrapper
from components.components_output_maker import CompleteFastaOutputMaker
from components.components_output_maker import CompleteFolderSummaryMaker
from components.components_output_maker import CompleteCasSummaryFolderMaker
from components.components_output_maker import CompleteJsonOutputMaker
from components.components_helpers import multiline_fasta_check, multiline_fasta_handle, multiline_fasta_handle_python
from components.components_helpers import folder_of_multifasta_handle

warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

FLAG_DEVELOPER_MODE = False

parser = argparse.ArgumentParser(description='Run Identifier')
parser.add_argument('--input_folder', type=str, default=None,
                    help='input folder (default: None)')

parser.add_argument('--file', type=str, default=None,
                    help='input file (default: None)')

parser.add_argument('--input_folder_multifasta', type=str, default=None,
                    help='input folder of multifasta (default: None)')

parser.add_argument('--model', type=str, default="ALL",
                    help='model_to_use (default: ALL)')

parser.add_argument('--additional_model', type=str, default=None,
                    help='model_to_use (default: None)')

parser.add_argument('--result_folder', type=str, default="Results",
                    help='folder with the result (default: Results)')

parser.add_argument('--pickle_report', type=str, default='',
                    help='pickled report file (default: None)')

parser.add_argument('--json_report', type=str, default='',
                    help='json report file (default: None)')

parser.add_argument('--fasta_report', type=str, default=False,
                    help='fasta report file (default: False)')

parser.add_argument('--strand', type=str, default=True,
                    help='CRISPR array orientation prediction (default: True)')

parser.add_argument('--cas', type=str, default=False,
                    help='cas genes computation (default: False)')

parser.add_argument('--is_element', type=str, default=True,
                    help='is element computation (default: True)')

parser.add_argument('--parallel', type=str, default=True,
                    help='parallel computations (default: True)')

parser.add_argument('--cpu', type=str, default="ALL",
                    help='parallel computations (default: ALL)')

parser.add_argument('--fast_run', type=str, default=False,
                    help='fast run option (default: False)')

parser.add_argument('--degenerated', type=bool, default=True,
                    help='degenerated_repeat_computation (default: True)')

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

parser.add_argument('--enhancement_max_min', type=bool, default=True,
                    help='enhancement with filter (default: True)')

parser.add_argument('--enhancement_start_end', type=bool, default=True,
                    help='enhancement with start end omitting (default: True)')

parser.add_argument('--max_identical_spacers', type=int, default=4,
                    help='maximum number of identical spacers in the array (default: 4)')

parser.add_argument('--max_identical_cluster_spacers', type=int, default=3,
                    help='maximum number of consecutive identical spacers in the array (default: 3)')

parser.add_argument('--margin_degenerated', type=int, default=30,
                    help='maximum length of the spacer margin for the degenerated search (default: 30)', )

parser.add_argument('--max_edit_distance_enhanced', type=int, default=6,
                    help='maximum edit distance for the evaluated array enhancement (default: 6)')


script_absolute_path = os.path.dirname(os.path.abspath(__file__))
work_directory = os.getcwd()
pid = os.getpid()

args = parser.parse_args()

complete_path_folder = (args.input_folder)
if complete_path_folder:
    complete_path_folder = Path(complete_path_folder).absolute()

complete_path_file = args.file
if complete_path_file:
    complete_path_file = Path(complete_path_file).absolute()

complete_folder_multifasta = args.input_folder_multifasta
if complete_folder_multifasta:
    complete_folder_multifasta = Path(complete_folder_multifasta).absolute()

folder_result = args.result_folder
if folder_result:
    folder_result = Path(folder_result).absolute()

pickle_folder = args.pickle_report
if pickle_folder:
    pickle_folder = Path(pickle_folder).absolute()

json_folder = args.json_report
if json_folder:
    json_folder = Path(json_folder).absolute()

list_models = ["8", "9", "10"] if args.model == "ALL" else [args.model]
flag_possible_differentiate_model = args.additional_model
if flag_possible_differentiate_model not in ["possible", "all"]:
    flag_possible_differentiate_model = None


flag_enhancement_max_min = args.enhancement_max_min
flag_enhancement_start_end = args.enhancement_start_end

flag_parallel = False if (args.parallel in ["False", False]) else True
flag_cpu = args.cpu
flag_fast_run = False if (args.fast_run in ["False", False]) else True

strand_flag = False if (args.strand in ["False", False]) else True
cas_flag = False if (args.cas in ["False", False]) else True
is_flag = False if (args.is_element in ["False", False]) else True
degenerated_flag = False if (args.degenerated in ["False", False]) else True
fasta_report = False if (args.fasta_report in ["False", False]) else True

flags = {"flag_parallel": flag_parallel,
         "flag_cpu": flag_cpu,
         "flag_fast_run": flag_fast_run,
         "flag_strand": strand_flag,
         "flag_cas": cas_flag,
         "flag_is": is_flag,
         "flag_fasta_report": fasta_report,
         "flag_degenerated": degenerated_flag,
         "flag_enhancement_min_max": flag_enhancement_max_min,
         "flag_enhancement_start_end": flag_enhancement_start_end
}

min_rep = args.min_len_rep
max_rep = args.max_len_rep
max_spacer = args.max_len_spacer
min_spacer = args.min_len_spacer
min_repeats = args.min_repeats
max_identical_spacers = args.max_identical_spacers
max_identical_cluster_spacers = args.max_identical_cluster_spacers
margin_degenerated = args.margin_degenerated
max_edit_distance_enhancement = args.max_edit_distance_enhanced

parameters = {
    "param_min_avg_repeat_length": min_rep,
    "param_max_avg_repeat_length": max_rep,
    "param_max_avg_spacer_length": max_spacer,
    "param_min_avg_spacer_length": min_spacer,
    "param_min_repeats": min_repeats,
    "param_max_identical_spacers": max_identical_spacers,
    "param_max_identical_cluster_spacers": max_identical_cluster_spacers,
    "param_spacer_margin_degenerated_search": margin_degenerated,
    "param_max_edit_distance": max_edit_distance_enhancement
}


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


pid_work_directory = os.path.join(work_directory, 'Identify_Temp' + str(pid))
if not os.path.exists(pid_work_directory):
    os.makedirs(pid_work_directory)
    os.chdir(pid_work_directory)


feature_list = ['.'.join([ALL_FEATURES[i] for i in best_combinations[model]]) for model in list_models]
list_ml_classifiers = [ClassifierWrapper(classifier_type=None,
                                         load_option=script_absolute_path + "/trained_models/extra_trees/extra_trees_subset{}features.pkl".
                                                     format(model))
                       for model in list_models]


def run_over_folder_of_files(folder, result_folder, pickle_folder, chunk_number=None, number_of_chunks=None):
    files = [f for f in listdir(folder) if isfile(join(folder, f))]
    files_name_fix = [f.replace("\r", "").replace("\t", "").replace("\n", "") for f in files]
    for old_name, new_name in zip(files, files_name_fix):
        old_path = join(folder, old_name)
        new_path = join(folder, new_name)
        if old_path != new_path:
            os.system(f"mv {old_path} {new_path}")
    files = sorted(files_name_fix)

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
        pl = Pipeline(result_folder_path="{}/".format(result_folder),
                      pickle_folder_path="{}".format(pickle_folder),
                      json_folder_path="{}".format(json_folder),
                      file_path=join(folder, file),
                      list_ml_classifiers=list_ml_classifiers,
                      list_features=feature_list,
                      parameters=parameters,
                      flags=flags,
                      flag_dev_mode=FLAG_DEVELOPER_MODE,
                      absolute_directory_path=script_absolute_path)

    cfsm = CompleteFolderSummaryMaker(folder_result=result_folder)
    ccfsm = CompleteCasSummaryFolderMaker(folder_result=result_folder)
    cfom = CompleteFastaOutputMaker(folder_result=result_folder)
    if json_folder:
        cjsm = CompleteJsonOutputMaker(folder_json_result=json_folder, folder_text_tesult=result_folder)


def run_over_one_file(file, result_folder, pickle_folder, json_folder):
    print("\n\n\n\t\t\t\tExecuting file {}\n\n\n".format(file))
    pl = Pipeline(result_folder_path="{}/".format(result_folder),
                  pickle_folder_path="{}".format(pickle_folder),
                  json_folder_path="{}".format(json_folder),
                  file_path=join(file),
                  list_ml_classifiers=list_ml_classifiers,
                  list_features=feature_list,
                  parameters=parameters,
                  flags=flags,
                  flag_dev_mode=FLAG_DEVELOPER_MODE,
                  absolute_directory_path=script_absolute_path)



def main():
    start_time = time()
    if complete_path_file:
        if multiline_fasta_check(complete_path_file):
            print("Multifasta")
            folder_multifasta = multiline_fasta_handle_python(complete_path_file)
            print(folder_multifasta)
            run_over_folder_of_files(folder_multifasta, folder_result, pickle_folder, json_folder)
            shutil.rmtree(folder_multifasta)
        else:
            run_over_one_file(complete_path_file, folder_result, pickle_folder, json_folder)
    elif complete_path_folder:
        run_over_folder_of_files(complete_path_folder, folder_result, pickle_folder, json_folder)
    elif complete_folder_multifasta:
        print("Folder Multifasta")
        folder_multifasta = folder_of_multifasta_handle(complete_folder_multifasta)
        run_over_folder_of_files(folder_multifasta, folder_result, pickle_folder, json_folder)
    else:
        print("No input was provided")

    end_time = time()
    print("Elapsed time: ", end_time-start_time)


if __name__ == "__main__":
    main()
    shutil.rmtree(pid_work_directory, ignore_errors=True)

