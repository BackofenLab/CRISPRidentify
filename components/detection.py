# coding: utf-8
"""The detection step of the pipeline"""

import os
from collections import namedtuple
from os import listdir
from os.path import isfile
from os.path import basename
from itertools import groupby
from functools import wraps
from multiprocessing import Pool
import multiprocessing

import regex
import json
from time import time

from crispr_candidate import CrisprCandidate
from components_filter_creation import FilterApproximation
from components_fuzzy_filter import AdvancedFuzzySearchFilter
from log_file_operations import logging

DEBUG_MODE = False
V_Repeat = namedtuple('V_Repeat', 'begin_first, begin_second, length, sequence')
ClusterSequence = namedtuple("ClusterSequence", "sequence, start, end, tuple_repeats")


class ClusterCandidate:
    def __init__(self, list_cluster_v_repeats, flag_fast_run=False,
                 flag_flag_enhancement_max_min=False, flag_enhancement_start_end=False):
        """
        Creates a class which has list of repeat candidates obtained from vmatch
        Repeat candidates are V_Repeat classes
        Has _calculate_different_repeat_sequences to calculate all different repeat candidate sequences
        """
        self.list_cluster_repeats = list_cluster_v_repeats
        self.begin = self.list_cluster_repeats[0].begin_first
        self.end = self.list_cluster_repeats[-1].begin_second + self.list_cluster_repeats[-1].length
        self.goes_through_zero = None

        self._calculate_different_repeat_sequences()
        if not flag_fast_run:
            if flag_flag_enhancement_max_min:
                self._enhance_with_filter_approximation()
            self._enhance_repeat_sequences()

            if flag_enhancement_start_end:
                self._adding_additional_candidates()

            self.list_clust_dif_rep_seq = sorted(self.list_clust_dif_rep_seq)

    def _calculate_different_repeat_sequences(self):
        self.list_clust_dif_rep_seq = list(set([r.sequence for r in self.list_cluster_repeats]))
        if DEBUG_MODE:
            print(self.begin)
            print(self.list_clust_dif_rep_seq)

    def _enhance_with_filter_approximation(self):
        fa = FilterApproximation(self.list_clust_dif_rep_seq)
        missing_repeats_from_filter = fa.output()

        self.list_clust_dif_rep_seq += missing_repeats_from_filter
        self.list_clust_dif_rep_seq = list(set(self.list_clust_dif_rep_seq))

    def _enhance_repeat_sequences(self):
        list_groups_of_repeats = []
        list_sorted_repeats_by_length = sorted(self.list_clust_dif_rep_seq, key=lambda x: len(x), reverse=True)
        for repeat_seq in list_sorted_repeats_by_length:
            if not list_groups_of_repeats:
                list_groups_of_repeats = [[repeat_seq]]
            else:
                self.add_new(list_groups_of_repeats, repeat_seq)

        list_complete_groups = [self.complete_with_all_intermediate(group_to_complete)
                                for group_to_complete in list_groups_of_repeats]

        list_all_repeats_candidates = [repeat_seq for group in list_complete_groups for repeat_seq in group]
        list_all_repeats_candidates = list(set(list_all_repeats_candidates))

        self.list_clust_dif_rep_seq = list_all_repeats_candidates
        if DEBUG_MODE:
            print(self.begin)
            print(self.list_clust_dif_rep_seq)

    def _adding_additional_candidates(self):
        copy_list_repeat_candidates = self.list_clust_dif_rep_seq[:]
        for repeat_seq in copy_list_repeat_candidates:
            for i in range(4):
                for j in range(4):
                    additional_candidate = repeat_seq[i:len(repeat_seq) - j]
                    if additional_candidate not in self.list_clust_dif_rep_seq:
                        self.list_clust_dif_rep_seq.append(additional_candidate)

        self.list_clust_dif_rep_seq = list(set(self.list_clust_dif_rep_seq))
        if DEBUG_MODE:
            print(self.begin)
            print(self.list_clust_dif_rep_seq)

    @staticmethod
    def add_new(groups, new_seq):
        for group_to_check in groups:
            if new_seq in group_to_check[0]:
                group_to_check.append(new_seq)
                break
        else:
            groups.append([new_seq])

    @staticmethod
    def complete_with_all_intermediate(group_to_complete):
        if len(group_to_complete) < 2:
            return group_to_complete
        else:
            group_to_complete = sorted(group_to_complete, key=lambda x: len(x), reverse=True)
            longest_one = group_to_complete[0]
            length_longest = len(longest_one)
            shortest_one = group_to_complete[-1]
            all_possible_substrings = [longest_one[i:j+1] for i in range(length_longest)
                                       for j in range(i, length_longest)
                                       if len(longest_one[i:j+1]) > len(shortest_one)]

            for candidate in all_possible_substrings:
                if shortest_one in candidate:
                    group_to_complete.append(candidate)

        return group_to_complete

    def __repr__(self):
        string = "Begin: {}\n".format(self.begin)
        string += "End: {}\n".format(self.end)
        string += "Repeat candidates:\n"
        for r in self.list_cluster_repeats:
            string += '{}\n'.format(r.sequence)

        return string


class FuzzySearch:
    def __init__(self, sequence, sequence_start, repeat, weighted_error):
        self.repeat_candidate = repeat
        self.pattern = "(?e)" + "({})".format(repeat) + weighted_error
        self.sequence = sequence
        self.sequence_start = sequence_start

        self._match()
        if self.match_hit is True:
            self._get_repeats()
            self._get_spacers()
            self._get_errors()
            self._get_absolute_starts()
            self._get_relative_error_indexes()
            self._get_all_insertions_deletions()
            self._get_gaped_repeats()
            self._compute_dot_representation()
        self.matching_object = None

    def _match(self):
        """Run fuzzy search on the sequence"""

        #self.matching_object = list(regex.finditer(self.pattern, self.sequence, ENHANCEMATCH=True))
        self.matching_object = list(regex.finditer(self.pattern, self.sequence))
        if self.matching_object:
            self.match_hit = True
        else:
            self.match_hit = False
        
    def _get_repeats(self):
        """Obtain the sequences that were matched"""
        self.list_repeats = [match.group() for match in self.matching_object]

    def _get_spacers(self):
        """Obtain the spacers sequences"""
        spacer_starts = [match.end() for match in self.matching_object][:-1]
        spacer_end = [match.start() for match in self.matching_object][1:]
        spacer_intervals = [(start, end) for start, end in zip(spacer_starts, spacer_end)]
        self.list_spacers = [self.sequence[interval[0]:interval[1]] for interval in spacer_intervals]

    def _get_errors(self):
        """Obtain all the errors"""
        self.number_substitutions = 0
        self.number_insertions = 0
        self.number_deletions = 0

        for match in self.matching_object:
            s, i, d = match.fuzzy_counts
            self.number_substitutions += s
            self.number_insertions += i
            self.number_deletions += d

        self.number_errors = self.number_substitutions + self.number_insertions + self.number_deletions

    def _get_absolute_starts(self):
        """Obtaining absolute start of matched sequences"""
        self.list_absolute_start = [self.sequence_start + match.start() for match in self.matching_object]
        self.start_end = '{}_{}'.format(self.list_absolute_start[0], self.list_absolute_start[-1])

    def _get_relative_error_indexes(self):
        """Obtaining indexes of errors which are related to each sequence and not global"""
        self.list_relative_error_indexes = []
        for match in self.matching_object:
            tuple_match_errors = match.fuzzy_changes
            list_relative_errors = [[e - match.start() for e in err_type] for err_type in tuple_match_errors]
            self.list_relative_error_indexes.append(list_relative_errors)

    def _get_all_insertions_deletions(self):
        self.list_all_insertions = []
        self.list_all_deletions = []

        for err_type in self.list_relative_error_indexes:
            self.list_all_insertions += self.unique_gaps_end(err_type[1])
            self.list_all_deletions += self.unique_gaps_end(err_type[2])

        self.list_all_insertions = list(set(self.list_all_insertions))
        self.list_all_deletions = list(set(self.list_all_deletions))

    @staticmethod
    def gap_introduction(sequence, list_of_gaps, gap_char):        
        new_seq = ""        
        for index in range(len(sequence) + len(list_of_gaps)):
            if index in list_of_gaps:
                new_seq += gap_char
            else:                
                new_seq += sequence[0]
                sequence = sequence[1:]
        
        return new_seq

    @staticmethod
    def apply_insertions_to_deletions(list_insertions, list_deletions):        
        list_deletions_new = []
        for element in list_deletions:
            add = sum(element >= i for i in list_insertions)
            list_deletions_new.append(element + add)
        return list_deletions_new
    
    @staticmethod
    def unique_gaps_end(list_gaps):
        list_gaps_new = []
        for el in list_gaps:
            while True:
                if el in list_gaps_new:
                    el += 1
                else:
                    list_gaps_new.append(el)
                    break
        return list_gaps_new
    
    @staticmethod
    def regex_bug_fix(list_gaps):
        return [x if x >= 0 else 0 for x in list_gaps]
    
    def _get_gaped_repeats(self):
        self.list_gaped_repeats = []        
        for index, repeat in enumerate(self.list_repeats):           
            local_insertions = self.list_relative_error_indexes[index][1]
            local_insertions = self.regex_bug_fix(local_insertions)
            local_insertions = self.unique_gaps_end(local_insertions)
            insertions = list(set(self.list_all_insertions).difference(set(local_insertions)))           
            
            deletions = self.list_relative_error_indexes[index][2]
            deletions = self.regex_bug_fix(deletions)
            deletions = self.unique_gaps_end(deletions)
            
            gaped_repeat = self.gap_introduction(repeat, deletions, "-")
            gaped_repeat = self.gap_introduction(gaped_repeat, insertions, " ")    
        
            self.list_gaped_repeats.append(gaped_repeat)
        self.gaped_matching_pattern = self.gap_introduction(self.repeat_candidate, self.list_all_insertions, " ")

    """
    def __repr__(self):
        string = ""
        if self.matching_object:
            string += repr(self.matching_object)
            string += '\n'
            for match in self.matching_object:
                string += str(match.fuzzy_changes)
            string += str(self.list_relative_error_indexes)
            string += '\n'
            max_length_start_index = max(len(str(start)) for start in self.list_absolute_start) + 3
            max_length_spacer = max(len(spacer) for spacer in self.list_spacers) + 3

            for index, gaped_repeat in enumerate(self.list_gaped_repeats):
                repeat_start_index = self.list_absolute_start[index] + 1
                n_gaps_after_start = max_length_start_index - len(str(self.list_absolute_start[index]))

                if index < len(self.list_spacers):
                    spacer = self.list_spacers[index]
                else:
                    spacer = ""
                n_gaps_after_spacer = max_length_spacer - len(spacer)

                s, i, d = self.matching_object[index].fuzzy_counts
                errors = "s:{} i:{} d{}".format(s, i, d)

                string += "{}{}{}  {}{}{}\n".format(repeat_start_index,
                                                    " " * n_gaps_after_start,
                                                    gaped_repeat,
                                                    spacer,
                                                    " " * n_gaps_after_spacer,
                                                    errors)

            string += "_" * 100 + "\n"       

            string += " " * max_length_start_index + self.gaped_matching_pattern
            string += " " * (max_length_spacer + 2) + "s:{} i:{} d:{}".format(self.number_substitutions,
                                                                              self.number_insertions,
                                                                              self.number_deletions) + "\n"

            string += "_" * 100 + "\n"
            string += " " * max_length_start_index + self.repeat_candidate      

        return string
        
    """

    def __repr__(self):
        return self.dot_representation
    
    @staticmethod
    def apply_dots(repeat_with_gaps, matching_pattern_with_gaps):
        dotted_repeat = ''
        if len(repeat_with_gaps) != len(matching_pattern_with_gaps):
            raise ValueError
        for char_repeat, char_pattern in zip(repeat_with_gaps, matching_pattern_with_gaps):
            if char_repeat == char_pattern:
                if char_repeat == ' ':
                    dotted_repeat += char_repeat
                else:
                    dotted_repeat += '.'
            else:
                dotted_repeat += char_repeat
        
        return dotted_repeat

    def _compute_dot_representation(self):
        string = ""
        if self.matching_object:
            max_length_start_index = max(len(str(start)) for start in self.list_absolute_start) + 3
            try:
                max_length_spacer = max(len(spacer) for spacer in self.list_spacers) + 3
            except ValueError:
                max_length_spacer = 0

            for index, gaped_repeat in enumerate(self.list_gaped_repeats):
                repeat_start_index = self.list_absolute_start[index] + 1
                n_gaps_after_start = max_length_start_index - len(str(self.list_absolute_start[index]))

                if index < len(self.list_spacers):
                    spacer = self.list_spacers[index]
                else:
                    spacer = ""
                n_gaps_after_spacer = max_length_spacer - len(spacer)

                s, i, d = self.matching_object[index].fuzzy_counts
                errors = "s:{} i:{} d{}".format(s, i, d)
                dotted_repeats = self.apply_dots(gaped_repeat, self.gaped_matching_pattern)

                string += "{}{}{}  {}{}{}{}\n".format(repeat_start_index,
                                                      " " * n_gaps_after_start,
                                                      dotted_repeats,
                                                      spacer,
                                                      " " * n_gaps_after_spacer,
                                                      errors,
                                                      self.list_relative_error_indexes[index])

            string += "_" * 100 + "\n"

            string += " " * max_length_start_index + self.gaped_matching_pattern
            string += " " * (max_length_spacer + 2) + "s:{} i:{} d:{}".format(self.number_substitutions,
                                                                              self.number_insertions,
                                                                              self.number_deletions) + "\n"

            string += "_" * 100 + "\n"
            string += " " * max_length_start_index + self.repeat_candidate

        self.dot_representation = string
    
    def dot_repr(self):
        return self.dot_representation


class Detector(object):
    def __init__(self, file_path, flag_parallel, flag_cpu, flag_fast_run, flag_enhancement_max_min,
                 flag_enhancement_start_end, parameters, log_file):
        self.file_path = file_path
        self.flag_parallel = flag_parallel
        self.flag_cpu = flag_cpu
        self.flag_fast_run = flag_fast_run
        self.flag_enhancement_max_min = flag_enhancement_max_min
        self.flag_enhancement_start_end = flag_enhancement_start_end
        self.param_min_avg_repeat_length = parameters["param_min_avg_repeat_length"]
        self.param_max_avg_repeat_length = parameters["param_max_avg_repeat_length"]
        self.param_min_avg_spacer_length = parameters["param_min_avg_spacer_length"]
        self.param_max_avg_spacer_length = parameters["param_max_avg_spacer_length"]
        self.param_min_repeats = parameters["param_min_repeats"]

        self.list_repeat_candidates = []
        self.list_repeat_candidates_old = []
        self.list_repeat_candidates_new = []
        self.list_clusters = []
        self.list_cluster_sequences = []
        self.settings = {}
        self.log_file = log_file

        if self.flag_fast_run:
            self._open_settings_fast_run()
        else:
            self._open_settings()

        self._create_folders()
        self._create_input_uppercase()
        self._get_double_dna()
        self._mkvtree_command()
        self._vmatch_command()
        self._take_vmatch_result_old()
        self._take_vmatch_results_new()
        self._check_old_new_vmatch()
        self._clean_after_mkv()
        self._find_clusters()

        if not self.flag_fast_run:
            self._merge_candidates_based_on_similarity()
        self._filter_cluster_duplicates()

        self._extract_sequences()

        if self.flag_parallel:
            self._fuzzy_search_parallel()
        else:
            self._fuzzy_search()

        self._write_fuzzy_candidates()

        self._filter_fuzzy_searches_same_start_end()
        self._apply_advanced_filters()
        
    def _open_settings(self):
        self.settings["vmatch -l"] = "21 -10 100"
        self.settings["vmatch -e"] = "1"
        self.settings["vmatch -evalue"] = "1"
        self.settings["flanking_region"] = "100"
        self.settings["max_length_between_repeats"] = "100"

    def _open_settings_fast_run(self):
        self.settings["vmatch -l"] = "23 25 60"
        self.settings["vmatch -e"] = "1"
        self.settings["vmatch -evalue"] = "1"
        self.settings["flanking_region"] = "100"
        self.settings["max_length_between_repeats"] = "150"

    def _create_folders(self):
        """ Method to create 2 folders:
        Result folder there the final result will be written.
        Temp folder is for intermediate computations.
        Creates three class attributes: file_base result_path and temp_path
        Takes: nothing
        Returns: nothing"""

        self.file_base = basename(self.file_path)

        """
        temp_path = "Temp/Temp_{}".format(self.file_base.split(".")[0])
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)
        
        
        self.temp_path = temp_path
        """

    def _create_input_uppercase(self):
        with open(self.file_path, 'r') as f:
            lines = f.readlines()

        self.input_header = lines[0]
        self.dna = ''.join([line.strip() for line in lines[1:]])

        with open("new_input.fa", "w") as f:
            f.write(self.input_header)
            upper_dna = self.dna.upper()
            f.write(upper_dna)

    def _get_double_dna(self):
        """Obtaining the header and the DNA sequence from the file
           Double DNA is obtained for detecting CRISPRs at the "end-start" region"""
        with open(self.file_path, 'r') as f:
            lines = f.readlines()

        self.input_header = lines[0]
        self.dna = ''.join([line.strip() for line in lines[1:]])
        self.dna_length = len(self.dna)
        self.dna += self.dna
        self.dna = self.dna.upper()
    
    def _mkvtree_command(self):
        """Executes the first necessary step for vmatch:
        mkvtree command
        Takes: nothing
        Returns: nothing"""

        print("Execute MKV TREE")
        cmd = "tools/vmatch/mkvtree -db new_input.fa -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1"
        os.system(cmd)
    
    def _vmatch_command(self):
        """Executes vmatch command with
        obtained from settings parameters"""

        print("Execute VMATCH")
        l_flag_val = self.settings["vmatch -l"]
        e_flag_val = self.settings["vmatch -e"]
        e_value_flag_val = self.settings["vmatch -evalue"]

        cmd = "tools/vmatch/vmatch " + "-l " + l_flag_val + " -evalue " + e_value_flag_val + " -e " + e_flag_val +\
              " -s leftseq " + " -absolute -nodist -noevalue -noscore -noidentity " + \
              "-sort ia -best 1000000 -selfun tools/vmatch/sel392.so 55 " + "new_input.fa" + " > " +\
              "vmatch_result.txt"

        os.system(cmd)

        cmd = "tools/vmatch/vmatch " + "-l " + l_flag_val + " -evalue " + e_value_flag_val + " -e " + e_flag_val + \
              " -s leftseq " + " -absolute -nodist -noevalue -noscore -noidentity " + \
              "-sort ia -best 1000000 " + "new_input.fa" + " > " + \
              "vmatch_result_new.txt"

        os.system(cmd)

    def _clean_after_mkv(self):
        """Cleans the folder moving the intermediate files into
        the corresponding temp folder"""

        onlyfiles = [f for f in listdir('.') if isfile(f)]
        to_move = []

        for el in onlyfiles:
            if "new_input.fa" in el:
                to_move.append(el)

        for el in to_move:
            os.remove(el)

        #os.remove("vmatch_result.txt")
        #os.remove("vmatch_result_new.txt")

    def _take_vmatch_result_old(self):
        """Works with the vmatch result file
        fills in self.list_repeat_candidates"""

        file_to_open = "vmatch_result.txt"
        list_param = []

        if self.log_file:
            logging(self.log_file, self.file_path+"\n")
            logging(self.log_file, "Vmatch Results\n")

        with open(file_to_open, "r") as f:
            for line in f:
                if line[0] == '>':
                    line = line.replace(">", "> ")
                    line = line.replace(' +', ' ')
                    list_param = line.split()
                elif len(line.rstrip()) > 0:
                    seq_dr = line.rstrip()
                    repeat = V_Repeat(int(list_param[2]), int(list_param[5]), int(list_param[1]), seq_dr)
                    self.list_repeat_candidates_old.append(repeat)

        if self.log_file:
            for repeat in self.list_repeat_candidates_old:
                line = f"{repeat.sequence}\t{repeat.begin_first}\t{repeat.begin_second}\t\t{repeat.length}\n"
                logging(self.log_file, line)
            logging(self.log_file, "\n")

    def _take_vmatch_results_new(self):
        cutoff = 55
        file_to_open = "vmatch_result_new.txt"

        if self.log_file:
            logging(self.log_file, self.file_path + "\n")
            logging(self.log_file, "Vmatch Results\n")

        with open(file_to_open, "r") as f:

            all_lines = f.readlines()
            indexes_headers = [index for index, line in enumerate(all_lines) if ">" in line]

            info_lines = all_lines[indexes_headers[0]:]

            list_param = None
            seq_dr = None

            for line in info_lines:
                if line[0] == '>':
                    if list_param:
                        if len(seq_dr) <= cutoff:
                            repeat = V_Repeat(int(list_param[2]), int(list_param[5]), int(list_param[1]), seq_dr)
                            self.list_repeat_candidates_new.append(repeat)
                            seq_dr = None

                    line = line.replace(">", "> ")
                    line = line.replace(' +', ' ')
                    list_param = line.split()

                elif len(line.rstrip()) > 0:
                    if not seq_dr:
                        seq_dr = line.rstrip()
                    else:
                        seq_dr += line.rstrip()

            if list_param:
                if len(seq_dr) <= cutoff:
                    repeat = V_Repeat(int(list_param[2]), int(list_param[5]), int(list_param[1]), seq_dr)
                    self.list_repeat_candidates_new.append(repeat)

        self.list_repeat_candidates = self.list_repeat_candidates_new

    def _check_old_new_vmatch(self):
        if len(self.list_repeat_candidates) != len(self.list_repeat_candidates_old):
            self.list_repeat_candidates = self.list_repeat_candidates_old
        else:
            if all([r_o == r_new for r_o, r_new in zip(self.list_repeat_candidates_old, self.list_repeat_candidates_new)]):
                self.list_repeat_candidates = self.list_repeat_candidates_new
            else:
                self.list_repeat_candidates = self.list_repeat_candidates_old

    def _find_clusters(self):
        """ Takes value of max allowed distance between repeats
        Separates found repeats into clusters accordingly
        by creating list of repeat candidates for each cluster.
        When distance between two repeats is larger than the threshold
        the cluster is written and creating of a new one is started."""

        print("Finding clusters")
        max_len_b_rep = int(self.settings["max_length_between_repeats"])
        list_repeat_candidates = []

        for repeat_cand in self.list_repeat_candidates:
            #print(repeat_cand)
            if not list_repeat_candidates:
                list_repeat_candidates.append(repeat_cand)
            elif repeat_cand.begin_first - list_repeat_candidates[-1].begin_first < max_len_b_rep:
                list_repeat_candidates.append(repeat_cand)
            else:
                self.list_clusters.append(ClusterCandidate(list_repeat_candidates,
                                                           self.flag_fast_run,
                                                           self.flag_enhancement_max_min,
                                                           self.flag_enhancement_start_end))
                list_repeat_candidates = [repeat_cand]

        if list_repeat_candidates:
            self.list_clusters.append(ClusterCandidate(list_repeat_candidates,
                                                       self.flag_fast_run,
                                                       self.flag_enhancement_max_min,
                                                       self.flag_enhancement_start_end))


        print("Clusters were found: got {} clusters".format(len(self.list_clusters)))

        if self.log_file:
            logging(self.log_file, "Original Clusters\n")
            for cluster in self.list_clusters:
                line = f"{cluster.begin}\t{cluster.end}\n"
                logging(self.log_file, line)
                for repeat in cluster.list_clust_dif_rep_seq:
                    logging(self.log_file, f"{repeat}\n")




    def _merge_candidates_based_on_similarity(self):
        cm = CandidateMerger(self.dna, self.list_clusters)
        self.list_clusters = cm.output()
        print("Clusters were merged: now got {} clusters".format(len(self.list_clusters)))
    
    def _filter_cluster_duplicates(self):
        for cluster_candidate in self.list_clusters:
            if cluster_candidate.begin < self.dna_length:
                if cluster_candidate.end > self.dna_length:
                    cluster_candidate.goes_through_zero = True
            else:
                self.list_clusters.remove(cluster_candidate)

        #print(self.list_clusters)

    def _extract_sequences(self):
        for cluster in self.list_clusters:
            seq_start = max(0, cluster.begin - int(self.settings["flanking_region"]))
            seq_end = min(len(self.dna), cluster.end + int(self.settings["flanking_region"]))
            cluster_seq = self.dna[seq_start:seq_end]
            tup_cluster_dif_rep = tuple(cluster.list_clust_dif_rep_seq)

            self.list_cluster_sequences.append(ClusterSequence(cluster_seq, seq_start, seq_end, tup_cluster_dif_rep))

    def _fuzzy_search(self):
        print("Running Fuzzy search")
        print("No parallel")
        if self.flag_fast_run:
            weighted_error = "{i<=3,d<=3,s<=3,i+d+s<=6}"
        else:
            weighted_error = "{i<=3,d<=3,s<=3,i+d+s<=6}"
        self.dict_search = {}

        for cluster_sequence in self.list_cluster_sequences:
            list_fuzzy_results = []
            for repeat in cluster_sequence.tuple_repeats:
                fuzzy_s = FuzzySearch(cluster_sequence.sequence, cluster_sequence.start,
                                      repeat, weighted_error)
                if fuzzy_s.match_hit:
                    if len(fuzzy_s.list_repeats) > 1:
                        list_fuzzy_results.append(fuzzy_s)

            self.dict_search[cluster_sequence] = list_fuzzy_results

    def _fuzzy_search_parallel(self):
        print("Running Fuzzy search")
        if self.flag_fast_run:
            weighted_error = "{i<=3,d<=3,s<=3,i+d+s<=6}"
        else:
            weighted_error = "{i<=3,d<=3,s<=3,i+d+s<=6}"
        self.dict_search = {}

        for cluster_sequence in self.list_cluster_sequences:
            nr = len(cluster_sequence.tuple_repeats)
            input_tuples = zip(cluster_sequence.tuple_repeats, [cluster_sequence.sequence] * nr,
                               [cluster_sequence.start] * nr, [weighted_error] * nr)

            num_workers_suggested = multiprocessing.cpu_count() if self.flag_cpu == "ALL" else int(self.flag_cpu)
            max_possible = multiprocessing.cpu_count()
            num_workers = num_workers_suggested if num_workers_suggested < max_possible else max_possible
            with Pool(num_workers) as p:
                fuzzy_results = p.map(self._parallel_run_fuzzy_run, input_tuples)
                fuzzy_results = [x for x in fuzzy_results if x.match_hit]
                fuzzy_results = [x for x in fuzzy_results if len(x.list_repeats) > 1]

            self.dict_search[cluster_sequence] = fuzzy_results

    @staticmethod
    def _parallel_run_fuzzy_run(input_tuple):
        repeat, sequence, start, weighted_error = input_tuple

        return FuzzySearch(sequence, start,
                           repeat, weighted_error)

    def _write_fuzzy_candidates(self):
        with open("all_candidates.out", "w") as f:
            for search_results in self.dict_search.values():
                for search in search_results:
                    f.write(repr(search))
                    f.write("\n\n")
                f.write("=" * 100)
                f.write("\n\n\n")
    
    def _filter_fuzzy_searches_same_start_end(self):
        """Here we filter fuzzy searches which have same start and end
           Meaning that they will produce the same consensus repeat and having different matching
           patterns does not make sense
           So we take only the best ones meaning that best - lowest number of errors"""

        self.dict_filtered_search_start_end = {}

        for cluster_seq, list_fuzzy_s in self.dict_search.items():
            list_start_end = [fuzzy_s.start_end for fuzzy_s in list_fuzzy_s]
            pattern_len = [len(fuzzy_s.repeat_candidate) for fuzzy_s in list_fuzzy_s]
            tuples_st_end_len = zip(list_start_end, pattern_len)

            list_categories = [[fuzzy_s for fuzzy_s in list_fuzzy_s if
                                (fuzzy_s.start_end, len(fuzzy_s.repeat_candidate)) == tuple_info]
                               for tuple_info in tuples_st_end_len]

            best_fuzzy_s = [sorted(category, key=lambda x: x.number_errors)[0]
                            for category in list_categories]

            best_fuzzy_s_unique_repeat = []
            u_repeats = []
            for b_fuz in best_fuzzy_s:
                repeat = b_fuz.repeat_candidate
                if repeat not in u_repeats:
                    u_repeats.append(repeat)
                    best_fuzzy_s_unique_repeat.append(b_fuz)
            
            self.dict_filtered_search_start_end[cluster_seq] = best_fuzzy_s_unique_repeat

    def _apply_advanced_filters(self):
        print("Filtering the candidates")
        afsf = AdvancedFuzzySearchFilter(min_column_dominance_repeat=0.4,
                                         max_spacer_length=140, max_column_dominance_spacer=0.8,
                                         max_allowed_conseq_spacers=3, max_allowed_same_spacers=4,
                                         max_inconsistent_columns=5,
                                         min_avg_repeat_length=self.param_min_avg_repeat_length,
                                         max_avg_repeat_length=self.param_max_avg_repeat_length,
                                         min_avg_spacer_length=self.param_min_avg_spacer_length,
                                         max_avg_spacer_length=self.param_max_avg_spacer_length,
                                         min_repeats=self.param_min_repeats)

        self.dict_advanced_filtered_fuzzy_search = {}
        for key, values in self.dict_filtered_search_start_end.items():
            list_filtered_advanced = [afsf(value) for value in values]
            list_filtered_advanced = [x for x in list_filtered_advanced if x]
            if not list_filtered_advanced:
                sorted_by_num_errors = sorted(list(values), key=lambda x: x.number_errors)
                if sorted_by_num_errors:
                    candidate_fewer_mismatches = sorted_by_num_errors[0]
                    self.dict_advanced_filtered_fuzzy_search[key] = [candidate_fewer_mismatches]
            else:
                self.dict_advanced_filtered_fuzzy_search[key] = list_filtered_advanced

    def output(self):
        dict_output = {}
        for key, list_fuzzy in self.dict_advanced_filtered_fuzzy_search.items():
            new_key = (key.start, key.end)
            list_crispr_candidates = [CrisprCandidate(fuzzy.list_repeats, fuzzy.list_gaped_repeats,
                                                      fuzzy.list_spacers, fuzzy.list_absolute_start)
                                      for fuzzy in list_fuzzy]
            
            dict_output[new_key] = list_crispr_candidates
            
        return dict_output


class CandidateMerger:
    def __init__(self, dna, list_clusters):
        self.dna = dna
        self.list_clusters = list_clusters

        #self._merge_by_distance()
        self._merge_base_on_similarity_to_repeat()

    def _merge_by_distance(self):
        list_merged_by_distance = []
        for cluster in self.list_clusters:
            if not list_merged_by_distance:
                list_merged_by_distance.append(cluster)
            else:
                last_cluster = list_merged_by_distance.pop()
                for el in self.merging_close(last_cluster, cluster):
                    list_merged_by_distance.append(el)

        self.list_clusters = list_merged_by_distance

    @staticmethod
    def merging_close(first_cluster, second_cluster):
        if abs(first_cluster.end - second_cluster.begin) <= 100:
            list_v_repeats_first = first_cluster.list_cluster_repeats
            list_cluster_second = second_cluster.list_cluster_repeats
            new_cluster_candidate = ClusterCandidate(list_v_repeats_first + list_cluster_second)
            return new_cluster_candidate,
        else:
            return first_cluster, second_cluster

    def _merge_base_on_similarity_to_repeat(self):
        flag_looping = True
        list_clusters = self.list_clusters
        while flag_looping:
            flag_looping, list_clusters = self.iterative_merging(self.dna, list_clusters)

        self.list_clusters = list_clusters

    @staticmethod
    def mergeable(dna, first_cluster, second_cluster):
        first_cluster_end = first_cluster.end
        second_cluster_begin = second_cluster.begin
        dna_interval = dna[int(first_cluster_end)-1: int(second_cluster_begin)+1]
        main_repeat_to_check = second_cluster.list_cluster_repeats[0].sequence
        main_repeat_to_check = main_repeat_to_check[3:-3]
        if len(dna_interval) < 500:
            if main_repeat_to_check:
                occurrences = dna_interval.count(main_repeat_to_check)
                if occurrences >= len(dna_interval) / 160:
                    return True
            if len(second_cluster.list_cluster_repeats) > 3:
                second_repeat_to_check = second_cluster.list_cluster_repeats[1].sequence
                second_repeat_to_check = second_repeat_to_check[3:-3]
                if second_repeat_to_check:
                    occurrences_second = dna_interval.count(second_repeat_to_check)
                    if occurrences_second >= len(dna_interval) / 160:
                        return True
                third_repeat_to_check = second_cluster.list_cluster_repeats[2].sequence
                third_repeat_to_check = third_repeat_to_check[3:-3]
                if third_repeat_to_check:
                    occurrences_third = dna_interval.count(third_repeat_to_check)
                    if occurrences_third >= len(dna_interval) / 180:
                        return True
        return False

    def iterative_merging(self, dna, list_clusters):
        temp_list_clusters = []
        for index, _ in enumerate(list_clusters[:-1]):
            first_cluster = list_clusters[index]
            second_cluster = list_clusters[index+1]
            if self.mergeable(dna, first_cluster, second_cluster):
                list_v_repeats_first = first_cluster.list_cluster_repeats
                list_cluster_second = second_cluster.list_cluster_repeats
                new_cluster_candidate = ClusterCandidate(list_v_repeats_first+list_cluster_second, flag_fast_run=False)
                temp_list_clusters.append(new_cluster_candidate)
                temp_list_clusters += list_clusters[index+2:]
                return True, temp_list_clusters
            else:
                temp_list_clusters.append(first_cluster)

        if list_clusters:
            temp_list_clusters.append(list_clusters[-1])
        return False, temp_list_clusters

    def output(self):
        return self.list_clusters



