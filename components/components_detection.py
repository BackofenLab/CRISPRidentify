import os
import subprocess
from collections import namedtuple
from os.path import basename
from os import listdir
from os.path import isfile
import regex


V_Repeat = namedtuple('V_Repeat', 'begin_first, begin_second, length, sequence')
ClusterSequence = namedtuple("ClusterSequence", "sequence, start, end, tuple_repeats")

#                    INITIAL REPEAT CANDIDATES
#########################################################################
#########################################################################


class VmatchRun:
    def __init__(self, file_path, flag_fast_run):
        self.file_path = file_path
        self.flag_fast_run = flag_fast_run
        self.list_repeat_candidates = []

        self._get_run_settings()
        self._create_folders()
        self._create_input_uppercase()
        self._mkvtree_command()
        self._vmatch_command()
        self._take_vmatch_results_new()
        self._clean_after_mkv()

    def _get_run_settings(self):
        if self.flag_fast_run:
            self.settings = {"vmatch -l": "23 25 60",
                             "vmatch -e": "1",
                             "vmatch -evalue": "1",
                             "flanking_region": "100",
                             "max_length_between_repeats": "100"}
        else:
            self.settings = {"vmatch -l": "21 -10 100",
                             "vmatch -e": "1",
                             "vmatch -evalue": "1",
                             "flanking_region": "100",
                             "max_length_between_repeats": "100"}

    def _create_folders(self):
        self.file_base = basename(self.file_path)

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
        with open(self.file_path, 'r') as f:
            lines = f.readlines()

        self.input_header = lines[0]
        self.dna = ''.join([line.strip() for line in lines[1:]])
        self.dna_length = len(self.dna)
        self.dna += self.dna
        self.dna = self.dna.upper()

    def _mkvtree_command(self):
        #cmd = "tools/vmatch/mkvtree -db new_input.fa -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1"
        no_binary_cmd = "mkvtree -db new_input.fa -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1"
        os.system(no_binary_cmd)

    def _vmatch_command(self):
        l_flag_val = self.settings["vmatch -l"]
        e_flag_val = self.settings["vmatch -e"]
        e_value_flag_val = self.settings["vmatch -evalue"]

        #cmd = "tools/vmatch/vmatch " + "-l " + l_flag_val + " -evalue " + e_value_flag_val + " -e " + e_flag_val + \
        #      " -s leftseq " + " -absolute -nodist -noevalue -noscore -noidentity " + \
        #      "-sort ia -best 1000000 " + "new_input.fa" + " > " + \
        #      "vmatch_result_new.txt"

        no_binary_cmd = "vmatch " + "-l " + l_flag_val + " -evalue " + e_value_flag_val + " -e " + e_flag_val + \
                         " -s leftseq " + " -absolute -nodist -noevalue -noscore -noidentity " + \
                         "-sort ia -best 1000000 " + "new_input.fa" + " > " + \
                         "vmatch_result_new.txt"

        os.system(no_binary_cmd)


    def _take_vmatch_results_new(self):
        cutoff = 55
        file_to_open = "vmatch_result_new.txt"
        with open(file_to_open, "r") as f:

            all_lines = f.readlines()
            indexes_headers = [index for index, line in enumerate(all_lines) if ">" in line]

            if indexes_headers:
                info_lines = all_lines[indexes_headers[0]:]
            else:
                info_lines = []

            for index, line in enumerate(info_lines):
                if line[0] == '>':
                    seq_dr = info_lines[index + 1].strip()
                    line = line.replace(">", "> ")
                    line = line.replace(' +', ' ')
                    list_param = line.split()
                    if len(seq_dr) >= cutoff:
                        seq_dr = seq_dr[:cutoff]
                        list_param[1] = cutoff

                    repeat = V_Repeat(int(list_param[2]), int(list_param[5]), int(list_param[1]), seq_dr)
                    self.list_repeat_candidates.append(repeat)

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
        try:
            os.remove("vmatch_result_new.txt")
        except Exception:
            pass

    def output(self):
        return self.list_repeat_candidates


#                       CLUSTERS
###########################################################################
###########################################################################


class ClusterCandidate:
    def __init__(self, list_cluster_v_repeats):

        self.list_cluster_repeats = list_cluster_v_repeats
        self.begin = self.list_cluster_repeats[0].begin_first
        self.end = self.list_cluster_repeats[-1].begin_second + self.list_cluster_repeats[-1].length
        self.goes_through_zero = None

        self.list_clust_dif_rep_seq = list(set([r.sequence for r in self.list_cluster_repeats]))

    def __repr__(self):
        string = f"Cluster {self.begin} - {self.end}: \n"
        for s in self.list_clust_dif_rep_seq:
            string += str(s) + "\n"
        return string


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
        dna_interval = dna[int(first_cluster_end) - 1: int(second_cluster_begin) + 1]
        if second_cluster_begin < first_cluster_end:
            return True
        elif len(dna_interval) < 100:
            candidate_first = set(first_cluster.list_clust_dif_rep_seq)
            candidate_second = set(second_cluster.list_clust_dif_rep_seq)
            overlap_candidates = candidate_first.intersection(candidate_second)
            if overlap_candidates:
                return True
        elif len(dna_interval) < 500:
            main_repeat_to_check = second_cluster.list_cluster_repeats[0].sequence
            main_repeat_to_check = main_repeat_to_check[3:-3]
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
                new_cluster_candidate = ClusterCandidate(list_v_repeats_first+list_cluster_second)
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


class ClusterMaker:
    def __init__(self, list_repeats_from_vmatch, dna):
        self.list_repeats_from_vmatch = list_repeats_from_vmatch
        self.dna = dna
        self.dna_length = len(dna)

        self.list_clusters = []
        self.list_cluster_sequences = []

        self._read_settings()
        self._find_clusters()
        self._merge_candidates_based_on_similarity()
        self._filter_cluster_duplicates()

    def _read_settings(self):
        self.settings = {"max_length_between_repeats": "100",
                         "flanking_region": "100"}

    def _find_clusters(self):
        max_len_b_rep = int(self.settings["max_length_between_repeats"])
        list_repeat_candidates = []

        for repeat_candidate in self.list_repeats_from_vmatch:
            if not list_repeat_candidates:
                list_repeat_candidates.append(repeat_candidate)
            elif repeat_candidate.begin_first - list_repeat_candidates[-1].begin_first < max_len_b_rep:
                list_repeat_candidates.append(repeat_candidate)
            else:
                self.list_clusters.append(ClusterCandidate(list_repeat_candidates))
                list_repeat_candidates = [repeat_candidate]

        if list_repeat_candidates:
            self.list_clusters.append(ClusterCandidate(list_repeat_candidates))

        #print("Clusters were found: got {} clusters".format(len(self.list_clusters)))

    def _merge_candidates_based_on_similarity(self):
        cm = CandidateMerger(self.dna, self.list_clusters)
        self.list_clusters = cm.output()
        #print("Clusters were merged: now got {} clusters".format(len(self.list_clusters)))

    def _filter_cluster_duplicates(self):
        for cluster_candidate in self.list_clusters:
            if cluster_candidate.begin < self.dna_length:
                if cluster_candidate.end > self.dna_length:
                    cluster_candidate.goes_through_zero = True
            else:
                self.list_clusters.remove(cluster_candidate)

    def output(self):
        return self.list_clusters


class FilterApproximationClusters:
    def __init__(self, clusters):
        self.clusters = clusters
        self.clusters_filter_enhanced = []

        self._apply_filter_enhancement()

    def _apply_filter_enhancement(self):
        for cluster in self.clusters:
            list_different_repeats = cluster.list_clust_dif_rep_seq
            fa = FilterApproximation(list_different_repeats)
            missing_repeats_from_filter = fa.output()

            list_different_repeats += missing_repeats_from_filter
            list_different_repeats = list(set(list_different_repeats))

            cluster.list_clust_dif_rep_seq = list_different_repeats
            self.clusters_filter_enhanced.append(cluster)

    def output(self):
        return self.clusters_filter_enhanced


class FilterApproximation:
    def __init__(self, list_repeats):
        self.list_repeats = list_repeats
        self.list_repeats = sorted(self.list_repeats)
        self.max_seq = None
        self.min_seq = None
        self.list_missing_candidates = []

        self._relative_path_generation()
        self._make_fasta_with_repeats()

        if len(self.list_repeats) > 2:
            self._obtain_max_sequence()
            self._obtain_min_sequence()
            self._compute_all_the_missing_cases()

        self._clean_up()

    def _relative_path_generation(self):
        full_path = os.path.realpath(__file__)
        self.absolute_path_to_tools = "/" + "/".join(full_path.split("/")[:-2]) + "/tools/"

    def _make_fasta_with_repeats(self):
        with open("clustal_repeats.fa", "w") as f:
            lines = "".join([f">r{index}\n{repeat}\n" for index, repeat in enumerate(self.list_repeats, 1)])
            f.write(lines)

    def _run_clustal_omega_repeats(self):
        #cmd = self.absolute_path_to_tools + "clustalOmega/clustalo -i clustal_repeats.fa -o loc_align.txt"
        cmd_no_binary = "clustalo -i clustal_repeats.fa -o loc_align.txt"
        process = subprocess.Popen(cmd_no_binary, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

    def _find_max_clustal_omega_repeat_results(self):
        with open("loc_align.txt") as f:
            self.align_result_repeats = f.readlines()
        if self.align_result_repeats:
            if set([line.strip()[0] for line in self.align_result_repeats[::2]]) == {">"}:
                alignments = [line.strip() for line in self.align_result_repeats[1::2]]
                columns = [[al[index] for al in alignments] for index in range(len(alignments[0]))]
                columns = [[x for x in column if x != "-"] for column in columns]
                columns = [sorted(column) for column in columns]
                # most_common = [max(set(lst).difference(set("-")), key=lst.count) for lst in columns]
                most_common = [max(lst, key=lst.count) for lst in columns]
                self.max_seq = "".join(most_common)
        os.remove("loc_align.txt")

    def _find_min_clustal_omega_repeat_results(self):
        index_start = None
        index_end = None

        with open("loc_align.txt") as f:
            self.align_result_repeats = f.readlines()
        if self.align_result_repeats:
            if set([line.strip()[0] for line in self.align_result_repeats[::2]]) == {">"}:
                alignments = [line.strip() for line in self.align_result_repeats[1::2]]
                columns = [[al[index] for al in alignments] for index in range(len(alignments[0]))]
                columns = [[x for x in column if x != "-"] for column in columns]
                columns = [sorted(column) for column in columns]
                most_common = [max(lst, key=lst.count) for lst in columns]
                conserved = "".join(most_common)

                for index in range(len(alignments[0])):
                    column = [al[index] for al in alignments]
                    if "-" not in column:
                        if not index_start:
                            index_start = index
                        index_end = index
                if index_start and index_end:
                    self.min_seq = conserved[index_start:index_end + 1]
        os.remove("loc_align.txt")

    def _run_clustal_omega_min_max(self):
        with open("min_max.fa", "w") as f:
            f.write(">max\n")
            f.write(f"{self.max_seq}\n")
            f.write(">min\n")
            f.write(f"{self.min_seq}\n")

        #cmd = self.absolute_path_to_tools + "clustalOmega/clustalo -i min_max.fa -o min_max_align.txt"
        cmd_no_binary = "clustalo -i min_max.fa -o min_max_align.txt"
        process = subprocess.Popen(cmd_no_binary, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

    def _obtain_max_sequence(self):
        self._run_clustal_omega_repeats()
        self._find_max_clustal_omega_repeat_results()

    def _obtain_min_sequence(self):
        self._run_clustal_omega_repeats()
        self._find_min_clustal_omega_repeat_results()

    def _compute_all_the_missing_cases(self):
        if self.max_seq and self.min_seq:
            if len(self.min_seq) > 4:
                if abs(len(self.max_seq) - len(self.min_seq)) < 16:
                    length_longest = len(self.max_seq)
                    length_shortest = len(self.min_seq)
                    all_possible_substrings = [self.max_seq[i:j + 1] for i in range(length_longest)
                                               for j in range(i, length_longest)
                                               if len(self.max_seq[i:j + 1]) > length_shortest]

                    for candidate in all_possible_substrings:
                        if self.min_seq in candidate:
                            self.list_missing_candidates.append(candidate)

    def _clean_up(self):
        try:
            os.remove("clustal_repeats.fa")
        except FileNotFoundError:
            pass

    def output(self):
        return self.list_missing_candidates


class StartEndEnhancementClusters:
    def __init__(self, clusters):
        self.clusters = clusters
        self.clusters_enhanced_start_end = []

        self._enhance_start_end()

    def _enhance_start_end(self):
        for cluster in self.clusters:
            additional_candidates = []
            list_different_repeats = cluster.list_clust_dif_rep_seq
            for repeat_seq in list_different_repeats:
                for i in range(4):
                    for j in range(4):
                        additional_candidate = repeat_seq[i:len(repeat_seq) - j]
                        additional_candidates.append(additional_candidate)

            enhanced_repeats = list(set(list_different_repeats + additional_candidates))
            cluster.list_clust_dif_rep_seq = enhanced_repeats
            self.clusters_enhanced_start_end.append(cluster)

    def output(self):
        return self.clusters_enhanced_start_end


class IntermediateEnhancementClusters:
    def __init__(self, clusters):
        self.clusters = clusters
        self.clusters_enhanced_intermediate = []

        self._enhance_intermediate()

    def _enhance_intermediate(self):
        for cluster in self.clusters:
            list_groups_of_repeats = []
            list_different_repeats = cluster.list_clust_dif_rep_seq
            list_sorted_repeats_by_length = sorted(list_different_repeats, key=lambda x: len(x), reverse=True)
            for repeat_seq in list_sorted_repeats_by_length:
                if not list_groups_of_repeats:
                    list_groups_of_repeats = [[repeat_seq]]
                else:
                    self.add_new(list_groups_of_repeats, repeat_seq)

            list_complete_groups = [self.complete_with_all_intermediate(group_to_complete)
                                    for group_to_complete in list_groups_of_repeats]

            list_all_repeats_candidates = [repeat_seq for group in list_complete_groups for repeat_seq in group]
            list_all_repeats_candidates = list(set(list_all_repeats_candidates))

            cluster.list_clust_dif_rep_seq = list_all_repeats_candidates
            self.clusters_enhanced_intermediate.append(cluster)

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

    def output(self):
        return self.clusters_enhanced_intermediate




#                           ARRAYS
###########################################################################
###########################################################################


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
        self.matching_object = list(regex.finditer(self.pattern, self.sequence))
        if self.matching_object:
            self.match_hit = True
        else:
            self.match_hit = False

    def _get_repeats(self):
        self.list_repeats = [match.group() for match in self.matching_object]

    def _get_spacers(self):
        spacer_starts = [match.end() for match in self.matching_object][:-1]
        spacer_end = [match.start() for match in self.matching_object][1:]
        spacer_intervals = [(start, end) for start, end in zip(spacer_starts, spacer_end)]
        self.list_spacers = [self.sequence[interval[0]:interval[1]] for interval in spacer_intervals]

    def _get_errors(self):
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
        self.list_absolute_start = [self.sequence_start + match.start() for match in self.matching_object]
        self.start_end = '{}_{}'.format(self.list_absolute_start[0], self.list_absolute_start[-1])

    def _get_relative_error_indexes(self):
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





