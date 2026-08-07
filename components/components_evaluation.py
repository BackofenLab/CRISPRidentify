import numpy as np
import math
import os
import re
import subprocess


from components.components_eden import load
from components.components_eden import fold as seq2graph
from components.components_eden import kernel_matrix


#        BULK FEATURE EXTRACTION
#########################################################
#########################################################

class BulkFeatureExtractorOrf:
    def __init__(self, dict_crispr_candidates):
        self.dict_crispr_candidates = dict_crispr_candidates

        self._create_repeat_intervals_for_ofr()
        self._extract_orf_scores()
        self._clean_up()

    def _create_repeat_intervals_for_ofr(self):
        self.dict_repeat_intervals_for_orf = {}
        for key, list_crispr_candidates in self.dict_crispr_candidates.items():
            for index, crispr_candidate in enumerate(list_crispr_candidates):
                list_repeats = crispr_candidate.list_repeats
                list_spacers = crispr_candidate.list_spacers

                list_length_repeats = [len(repeat) for repeat in list_repeats]
                list_length_spacers = [len(spacer) for spacer in list_spacers]

                list_sum_list_repeats = [sum(list_length_repeats[:x]) for x in range(len(list_length_repeats) + 1)]
                list_sum_list_spacers = [sum(list_length_spacers[:x]) for x in range(len(list_length_spacers) + 1)]

                list_repeat_starts = [r + s for r, s in zip(list_sum_list_repeats, list_sum_list_spacers)]
                list_repeat_ends = [r + s for r, s in zip(list_sum_list_repeats[1:], list_sum_list_spacers)]

                list_repeat_starts = [r_start + 1 for r_start in list_repeat_starts]
                list_repeat_ends = [r_end + 1 for r_end in list_repeat_ends]

                list_repeat_intervals = [(r_start, r_end) for
                                         r_start, r_end in zip(list_repeat_starts, list_repeat_ends)]

                if key in self.dict_repeat_intervals_for_orf:
                    self.dict_repeat_intervals_for_orf[key].append(list_repeat_intervals)
                else:
                    self.dict_repeat_intervals_for_orf[key] = [list_repeat_intervals]

    def _extract_orf_scores(self):
        with open("file_for_prodigal.fa", "w") as f:
            for key, list_crispr_candidates in self.dict_crispr_candidates.items():
                for index, crispr_candidate in enumerate(list_crispr_candidates):
                    list_repeats = crispr_candidate.list_repeats
                    list_spacers = crispr_candidate.list_spacers

                    crispr_seq = ""
                    for r, s in zip(list_repeats, list_spacers):
                        crispr_seq += (r + s)
                    crispr_seq += list_repeats[-1]

                    f.write(f'>{key}_{index}\n')
                    f.write(crispr_seq)
                    f.write("\n")

        #cmd = 'tools/prodigal/prodigal'
        #cmd += ' -i file_for_prodigal.fa -o prodigal_result.txt -c -p meta'

        no_binary_cmd = "prodigal -i file_for_prodigal.fa -o prodigal_result.txt -c -p meta"

        process = subprocess.Popen(no_binary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        with open("prodigal_result.txt") as f:
            self.dict_orf_scores_different_inputs = {}
            lines_prodigal = f.readlines()

            for index, line in enumerate(lines_prodigal):
                if 'CDS' in line:
                    header_line = ""
                    for i in range(1, 100):
                        possible_header = lines_prodigal[index - i]
                        if "DEFINITION" in possible_header:
                            header_line = possible_header
                            break

                    if header_line:
                        interval_index_tuple = header_line.split("seqhdr=\"")[1].split("\";")[0]

                        next_line = lines_prodigal[index + 1]
                        region_string = re.search(r'\d+\.\.\d+', line).group(0)
                        region_tuple = tuple(int(x) for x in region_string.split('..'))
                        flag_complement = 1 if 'complement' in line else 0
                        tuple_region_flag_complement = region_tuple[0], region_tuple[1], flag_complement

                        match = re.search(r'conf=\d+\.\d+;', next_line)
                        confidence = float(match.group(0)[5:-1])

                        if interval_index_tuple in self.dict_orf_scores_different_inputs:
                            self.dict_orf_scores_different_inputs[interval_index_tuple][tuple_region_flag_complement] = confidence
                        else:
                            self.dict_orf_scores_different_inputs[interval_index_tuple] = {tuple_region_flag_complement: confidence}

        self.dict_final_orf_result = {}

        for key, list_crispr_candidates in self.dict_crispr_candidates.items():
            list_scores = []
            for index, crispr_candidate in enumerate(list_crispr_candidates):
                interval_index_tuple = f"{key}_{index}"
                if interval_index_tuple in self.dict_orf_scores_different_inputs:
                    orf_best_score = 0.0

                    list_repeats = crispr_candidate.list_repeats
                    list_spacers = crispr_candidate.list_spacers
                    crispr_length = sum(len(x) for x in (list_repeats + list_spacers))

                    list_tuple_region_flag_complement = self.dict_orf_scores_different_inputs[interval_index_tuple].keys()
                    for tuple_region_flag_complement in list_tuple_region_flag_complement:
                        if tuple_region_flag_complement[2] == 0:
                            reg_start = tuple_region_flag_complement[0]
                        else:
                            reg_start = crispr_length - tuple_region_flag_complement[1]

                        list_repeat_intervals = self.dict_repeat_intervals_for_orf[key][index]
                        for repeat_interval in list_repeat_intervals:
                            interval_start = repeat_interval[0]
                            interval_end = repeat_interval[1]

                            if interval_end >= reg_start >= interval_start:
                                score = self.dict_orf_scores_different_inputs[interval_index_tuple][tuple_region_flag_complement]
                                if score > orf_best_score:
                                    orf_best_score = score

                    list_scores.append(orf_best_score)

                else:
                    list_scores.append(0.0)

            self.dict_final_orf_result[key] = list_scores

    def _clean_up(self):
        try:
            os.remove("file_for_prodigal.fa")
        except OSError:
            pass

        try:
            pass
            os.remove("prodigal_result.txt")
        except OSError:
            pass

    def output_ofr_scores(self):
        return self.dict_final_orf_result


class BulkFeatureExtractorORFHMMR:
    def __init__(self, dict_crispr_candidates):
        self.dict_crispr_candidates = dict_crispr_candidates

        self._create_repeat_intervals_for_ofr()
        self._extract_orf_scores_and_proteins()
        self._run_hmm_search()
        self._extract_best_score_from_hmm()
        self._clean_up()

    def _create_repeat_intervals_for_ofr(self):
        self.dict_repeat_intervals_for_orf = {}
        for key, list_crispr_candidates in self.dict_crispr_candidates.items():
            for index, crispr_candidate in enumerate(list_crispr_candidates):
                list_repeats = crispr_candidate.list_repeats
                list_spacers = crispr_candidate.list_spacers

                list_length_repeats = [len(repeat) for repeat in list_repeats]
                list_length_spacers = [len(spacer) for spacer in list_spacers]

                list_sum_list_repeats = [sum(list_length_repeats[:x]) for x in range(len(list_length_repeats) + 1)]
                list_sum_list_spacers = [sum(list_length_spacers[:x]) for x in range(len(list_length_spacers) + 1)]

                list_repeat_starts = [r + s for r, s in zip(list_sum_list_repeats, list_sum_list_spacers)]
                list_repeat_ends = [r + s for r, s in zip(list_sum_list_repeats[1:], list_sum_list_spacers)]

                list_repeat_starts = [r_start + 1 for r_start in list_repeat_starts]
                list_repeat_ends = [r_end + 1 for r_end in list_repeat_ends]

                list_repeat_intervals = [(r_start, r_end) for
                                         r_start, r_end in zip(list_repeat_starts, list_repeat_ends)]

                if key in self.dict_repeat_intervals_for_orf:
                    self.dict_repeat_intervals_for_orf[key].append(list_repeat_intervals)
                else:
                    self.dict_repeat_intervals_for_orf[key] = [list_repeat_intervals]

    def _extract_orf_scores_and_proteins(self):
        with open("file_for_prodigal.fa", "w") as f:
            for key, list_crispr_candidates in self.dict_crispr_candidates.items():
                for index, crispr_candidate in enumerate(list_crispr_candidates):
                    list_repeats = crispr_candidate.list_repeats
                    list_spacers = crispr_candidate.list_spacers

                    crispr_seq = ""
                    for r, s in zip(list_repeats, list_spacers):
                        crispr_seq += (r + s)
                    crispr_seq += list_repeats[-1]

                    start, end = key
                    f.write(f">{start}_{end}_{index}\n")
                    f.write(crispr_seq)
                    f.write("\n")

        #cmd = 'tools/prodigal/prodigal'
        #cmd += ' -i file_for_prodigal.fa -o prodigal_result.txt -c -p meta'

        no_binary_cmd = "prodigal -i file_for_prodigal.fa -o prodigal_result.txt -c -p meta"
        process = subprocess.Popen(no_binary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        #cmd = 'tools/prodigal/prodigal'
        #cmd += ' -i file_for_prodigal.fa -p meta -a protein_results.fa'

        no_binary_cmd = "prodigal  -i file_for_prodigal.fa -p meta -a protein_results.fa"

        process = subprocess.Popen(no_binary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        with open("prodigal_result.txt") as f:
            self.dict_orf_scores_different_inputs = {}
            lines_prodigal = f.readlines()

            for index, line in enumerate(lines_prodigal):
                if 'CDS' in line:
                    header_line = ""
                    for i in range(1, 100):
                        possible_header = lines_prodigal[index - i]
                        if "DEFINITION" in possible_header:
                            header_line = possible_header
                            break

                    if header_line:
                        interval_index_tuple = header_line.split("seqhdr=\"")[1].split("\";")[0]

                        next_line = lines_prodigal[index + 1]
                        region_string = re.search(r'\d+\.\.\d+', line).group(0)
                        region_tuple = tuple(int(x) for x in region_string.split('..'))
                        flag_complement = 1 if 'complement' in line else 0
                        tuple_region_flag_complement = region_tuple[0], region_tuple[1], flag_complement

                        match = re.search(r'conf=\d+\.\d+;', next_line)
                        confidence = float(match.group(0)[5:-1])

                        if interval_index_tuple in self.dict_orf_scores_different_inputs:
                            self.dict_orf_scores_different_inputs[interval_index_tuple][
                                tuple_region_flag_complement] = confidence
                        else:
                            self.dict_orf_scores_different_inputs[interval_index_tuple] = {
                                tuple_region_flag_complement: confidence}

        self.dict_final_orf_result = {}

        for key, list_crispr_candidates in self.dict_crispr_candidates.items():
            list_scores = []
            for index, crispr_candidate in enumerate(list_crispr_candidates):
                start, end = key
                interval_index_tuple = f"{start}_{end}_{index}"
                if interval_index_tuple in self.dict_orf_scores_different_inputs:
                    orf_best_score = 0.0

                    list_repeats = crispr_candidate.list_repeats
                    list_spacers = crispr_candidate.list_spacers
                    crispr_length = sum(len(x) for x in (list_repeats + list_spacers))

                    list_tuple_region_flag_complement = self.dict_orf_scores_different_inputs[
                        interval_index_tuple].keys()
                    for tuple_region_flag_complement in list_tuple_region_flag_complement:
                        if tuple_region_flag_complement[2] == 0:
                            reg_start = tuple_region_flag_complement[0]
                        else:
                            reg_start = crispr_length - tuple_region_flag_complement[1]

                        list_repeat_intervals = self.dict_repeat_intervals_for_orf[key][index]
                        for repeat_interval in list_repeat_intervals:
                            interval_start = repeat_interval[0]
                            interval_end = repeat_interval[1]

                            if interval_end >= reg_start >= interval_start:
                                score = self.dict_orf_scores_different_inputs[interval_index_tuple][
                                    tuple_region_flag_complement]
                                if score > orf_best_score:
                                    orf_best_score = score

                    list_scores.append(orf_best_score)

                else:
                    list_scores.append(0.0)

            self.dict_final_orf_result[key] = list_scores

    def _run_hmm_search(self):
        #cmd = 'tools/hmm_search/hmmsearch --tblout result_hmm.out tools/hmm_search/models_tandem.hmm protein_results.fa'

        no_binary_cmd = "hmmsearch --tblout result_hmm.out ../tools/hmm_search/models_tandem.hmm protein_results.fa"

        process = subprocess.Popen(no_binary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

    def _extract_best_score_from_hmm(self):
        self.dict_hmm_results = {}
        self.dict_final_hmm_results = {}

        for key, list_crispr_candidates in self.dict_crispr_candidates.items():
            self.dict_hmm_results[key] = [None] * len(list_crispr_candidates)

        try:
            with open("result_hmm.out") as f:
                lines = f.readlines()

                for line in lines[3:-9]:
                    list_info = line.split()
                    if len(list_info) > 6:
                        score = float(list_info[5])
                        start, end, index = list_info[0].split("_")[:3]
                        index = int(index)

                        key = int(start), int(end)
                        if self.dict_hmm_results[key][index] is None:
                            self.dict_hmm_results[key][index] = score
                        else:
                            old_score = self.dict_hmm_results[key][index]
                            if score > old_score:
                                self.dict_hmm_results[key][index] = score

        except FileNotFoundError:
            pass

        for list_scores in self.dict_hmm_results.values():
            for index, value in enumerate(list_scores):
                if value is None:
                    list_scores[index] = 0.0

        for key in self.dict_crispr_candidates.keys():
            self.dict_final_hmm_results[key] = self.dict_hmm_results[key]

    def _clean_up(self):
        try:
            pass
            os.remove("file_for_prodigal.fa")
        except OSError:
            pass

        try:
            pass
            os.remove("prodigal_result.txt")
        except OSError:
            pass

        try:
            pass
            os.remove("protein_results.fa")
        except OSError:
            pass

        try:
            pass
            os.remove("result_hmm.out")
        except OSError:
            pass

    def output_ofr_scores(self):
        return self.dict_final_orf_result

    def output_hmmr_scores(self):
        return self.dict_final_hmm_results


class BulkFeatureExtractorBlast:
    def __init__(self, dict_crispr_candidates):
        self.dict_crispr_candidates = dict_crispr_candidates

        self.dict_blast_scores_1 = {}
        self.dict_blast_scores_2 = {}

        self._extract_all_blast_scores()
        self._clean_up()

    def _extract_all_blast_scores(self):
        list_consensus_rep = []
        for key, list_crispr_candidates in self.dict_crispr_candidates.items():
            for index, crispr_candidate in enumerate(list_crispr_candidates):
                list_consensus_rep.append(crispr_candidate.consensus)

        with open('file_with_all_consensus.fa', 'w') as f:
            for index, consensus_seq in enumerate(list_consensus_rep, 1):
                f.write(f">{index}\n")
                f.write(f"{consensus_seq}\n")

        db_file = 'Verified_repeats_dataset1.fa'

        #cmd = 'tools/blasting/blastn -query file_with_all_consensus.fa'
        #cmd += ' -db tools/blasting/'
        #cmd += db_file
        #cmd += ' -word_size 6'
        #cmd += ' -outfmt 6 -out output_fasta_bulk_extraction1'

        no_binary_cmd = f"blastn -query file_with_all_consensus.fa -db ../tools/blasting/{db_file} -word_size 6  -outfmt 6 -out output_fasta_bulk_extraction1"

        process = subprocess.Popen(no_binary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        db_file = 'Verified_repeats_dataset2.fa'

        #cmd = 'tools/blasting/blastn -query file_with_all_consensus.fa'
        #cmd += ' -db tools/blasting/'
        #cmd += db_file
        #cmd += ' -word_size 6'
        #cmd += ' -outfmt 6 -out output_fasta_bulk_extraction2'

        no_binary_cmd = f"blastn -query file_with_all_consensus.fa -db ../tools/blasting/{db_file} -word_size 6  -outfmt 6 -out output_fasta_bulk_extraction2"

        process = subprocess.Popen(no_binary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        with open("output_fasta_bulk_extraction1") as f:
            lines = f.readlines()
            max_num = len(list_consensus_rep)
            first_lines_scores = {}

            for line in lines:
                input_index = int(line.split()[0])
                if input_index not in first_lines_scores:
                    first_lines_scores[input_index] = float(line.strip().split()[-1])

            scores1 = []

            for possible_input_input_index in range(1, max_num+1):
                if possible_input_input_index in first_lines_scores:
                    bit_score = first_lines_scores[possible_input_input_index]
                else:
                    bit_score = 0.0

                scores1.append(bit_score)

        with open("output_fasta_bulk_extraction2") as f:
            lines = f.readlines()
            max_num = len(list_consensus_rep)
            first_lines_scores = {}

            for line in lines:
                input_index = int(line.split()[0])
                if input_index not in first_lines_scores:
                    first_lines_scores[input_index] = float(line.strip().split()[-1])

            scores2 = []

            for possible_input_input_index in range(1, max_num + 1):
                if possible_input_input_index in first_lines_scores:
                    bit_score = first_lines_scores[possible_input_input_index]
                else:
                    bit_score = 0.0

                scores2.append(bit_score)

        start_number = 0
        for key, list_crispr_candidates in self.dict_crispr_candidates.items():
            length_candidates = len(list_crispr_candidates)
            end_number = start_number + length_candidates
            scores_to_take1 = scores1[start_number:end_number]
            scores_to_take2 = scores2[start_number:end_number]
            self.dict_blast_scores_1[key] = scores_to_take1
            self.dict_blast_scores_2[key] = scores_to_take2
            start_number = end_number

    def _clean_up(self):
        try:
            os.remove("file_with_all_consensus.fa")
        except OSError:
            pass

        try:
            os.remove("output_fasta_bulk_extraction1")
        except OSError:
            pass

        try:
            os.remove("output_fasta_bulk_extraction2")
        except OSError:
            pass

    def output_blast(self):
        return self.dict_blast_scores_1, self.dict_blast_scores_2


class BulkFeatureExtractorMFE:
    def __init__(self, dict_crispr_candidates):
        self.dict_crispr_candidates = dict_crispr_candidates
        self.dict_final_mfe_result = {}

        self._extract_mfe_scores()
        self._clean_up()

    def _extract_mfe_scores(self):

        with open("file_for_mfe.fa", "w") as f:
            for key, list_crispr_candidates in self.dict_crispr_candidates.items():
                for index, crispr_candidate in enumerate(list_crispr_candidates):
                    consensus = crispr_candidate.consensus

                    f.write(f'>{key}_{index}\n')
                    f.write(consensus)
                    f.write("\n")

        #cmd = "cat file_for_mfe.fa | tools/rna_fold/RNAfold --noLP --noPS > rna_fold_output.txt"

        no_binary_cmd = "cat file_for_mfe.fa | RNAfold --noLP --noPS > rna_fold_output.txt"

        process = subprocess.Popen(no_binary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        with open("rna_fold_output.txt") as f:
            lines = f.readlines()

        header_lines = lines[::3]
        headers = [line.strip()[1:] for line in header_lines]
        score_lines = lines[2::3]
        scores = []
        for score_line in score_lines:
            try:
                search = re.search(r'-*\d+.\d+', score_line)
                mfe_score = search.group(0)
            except IndexError:
                mfe_score = 0.0

            scores.append(mfe_score)

        for header, score in zip(headers, scores):
            interval = header.split("_")[0]
            interval = interval[1:-1]
            interval = int(interval.split(", ")[0]), int(interval.split(", ")[1])
            if interval in self.dict_final_mfe_result:
                self.dict_final_mfe_result[interval].append(score)
            else:
                self.dict_final_mfe_result[interval] = [score]

    def _clean_up(self):
        try:
            os.remove("file_for_mfe.fa")
        except OSError:
            pass

        try:
            os.remove("rna_fold_output.txt")
        except OSError:
            pass

    def output_mfe_scores(self):
        return self.dict_final_mfe_result


class BulkFeatureExtractor:
    def __init__(self, dict_crispr_candidates):
        bfe_blast = BulkFeatureExtractorBlast(dict_crispr_candidates)
        bfe_orf_hmmr_score = BulkFeatureExtractorORFHMMR(dict_crispr_candidates)
        bfe_mfe = BulkFeatureExtractorMFE(dict_crispr_candidates)

        self.results_blast_scores = bfe_blast.output_blast()
        self.results_orf_scores = bfe_orf_hmmr_score.output_ofr_scores()
        self.results_hmm_results = bfe_orf_hmmr_score.output_hmmr_scores()
        self.result_mfe_scores = bfe_mfe.output_mfe_scores()

    def output(self):
        return self.results_blast_scores, self.results_orf_scores,\
               self.results_hmm_results, self.result_mfe_scores


#            REGULAR FEATURE EXTRACTION
#####################################################
#####################################################
class ConsensusFeatures:
    def __init__(self, index, consensus_seq):
        self.index = index
        self.consensus_seq = consensus_seq

    def output(self):
        raise NotImplementedError


class RepeatFeatures:
    def __init__(self, index, list_repeats):
        self.index = index
        self.list_repeats = list_repeats

    def output(self):
        raise NotImplementedError


class SpacerFeatures:
    def __init__(self, index, list_spacers):
        self.index = index
        self.list_spacers = list_spacers

    def output(self):
        raise NotImplementedError


class RepeatSpacersFeatures:
    def __init__(self, index, list_repeats, list_spacers):
        self.index = index
        self.list_repeats = list_repeats
        self.list_spacers = list_spacers

    def output(self):
        raise NotImplementedError


class CrisprNumberRepeats(RepeatFeatures):
    def output(self):
        number_repeats = len(self.list_repeats)
        return float(number_repeats)


class CrisprRepeatLen(RepeatFeatures):
    def output(self):
        return float(max((len(repeat) for repeat in self.list_repeats)))


class CrisprNumberMismatches(RepeatFeatures):
    def __init__(self, index, list_repeats, consensus_repeat_gaped):
        super().__init__(index, list_repeats)
        self.num_mismatches = sum(c_rep != c_c_repeat for repeat in list_repeats
                                  for c_rep, c_c_repeat in zip(repeat, consensus_repeat_gaped))

    def output(self):
        return float(self.num_mismatches)


class CrisprAvgSpacerLength(SpacerFeatures):
    def output(self):
        self.list_spacers = self.list_spacers
        if self.list_spacers[-1] == '':
            self.list_spacers = self.list_spacers[:-1]
        try:
            avg_spacer_length = sum(len(spacer) for spacer in self.list_spacers) / float(len(self.list_spacers))
        except ZeroDivisionError:
            avg_spacer_length = 0.0
        return float(avg_spacer_length)


class CrisprATRich(RepeatSpacersFeatures):
    def __init__(self, index, list_repeats, list_spacers):
        super().__init__(index, list_repeats, list_spacers)

        self.at_richness = None
        self._compute_at_rich()

    def _compute_at_rich(self):
        at_number = 0
        gc_number = 0
        for repeat in self.list_repeats:
            for char in repeat:
                if char == 'A' or char == 'T':
                    at_number += 1
                elif char == 'G' or char == 'C':
                    gc_number += 1

        try:
            self.at_richness = float(at_number) / (at_number + gc_number)
        except ZeroDivisionError:
            self.at_richness = 0.0

    def output(self):
        return float(self.at_richness)


class CrisprSpacerEveness(SpacerFeatures):
    def __init__(self, index, list_spacers, list_right_boundaries=(15, 25, 36, 52, 72)):

        super().__init__(index, list_spacers)

        self.list_right_boundaries = list_right_boundaries

        self.spacer_eveness = None
        self._compute_spacer_eveness()

    def _compute_spacer_eveness(self):
        list_bin_values = [0 for _ in range((len(self.list_right_boundaries) + 1))]
        for spacer in self.list_spacers:
            for boundary_index, right_boundary in enumerate(self.list_right_boundaries):
                if len(spacer) <= right_boundary:
                    list_bin_values[boundary_index] += 1
                    break
            else:
                list_bin_values[-1] += 1

        list_bin_normalized_values = [x / float(sum(list_bin_values)) for x in list_bin_values]
        summ = sum([(-1) * x * math.log(x + 0.00000001) for x in list_bin_normalized_values])
        self.spacer_evenness = round(summ / (math.log(len(list_bin_values))), 6)

    def output(self):
        return float(self.spacer_evenness)


class CrisprSimilarity(RepeatSpacersFeatures):
    def __init__(self, index, list_repeats, list_spacers):

        super().__init__(index, list_repeats, list_spacers)

        if self.list_spacers:
            if self.list_spacers[-1] == '':
                self.list_spacers = self.list_spacers[:-1]

        if set(self.list_spacers) == {"-"}:
            self.list_spacers = []

        self.similarity_score_repeats = None
        self.similarity_score_spacers = None
        self._compute_similarity_repeats_spacers()

    @staticmethod
    def check_list(input_list):
        for element in input_list:
            if len(set(element) - {"", "-"}) > 0:
                return True
        return False

    def _write_input_file_for_similarity(self, list_of_sequences):
        with open('file_for_similarity_{}.fa'.format(self.index), 'w') as f:
            for element_index, element in enumerate(list_of_sequences):
                f.write('>{}\n'.format(element_index + 1))
                f.write('{}\n'.format(element))

    def _make_data(self, fname):
        seqs = list(load(fname))
        pos_graphs = list(seq2graph(seqs))
        graphs = pos_graphs
        targets = np.array([1] * len(pos_graphs))

        return graphs, targets

    def _compute_similarity(self):
        graphs, targets = self._make_data('file_for_similarity_{}.fa'.format(self.index))
        os.remove('file_for_similarity_{}.fa'.format(self.index))

        x = kernel_matrix(graphs, r=3, d=4)
        n = x.shape[0]

        if n == 1:
            return 1.0
        else:
            return (np.sum(x) - n) / (n ** 2 - n)

    def _compute_similarity_repeats_spacers(self):
        if self.list_repeats:
            if self.check_list(self.list_repeats):
                self._write_input_file_for_similarity(self.list_repeats)
                self.similarity_score_repeats = self._compute_similarity()
                if self.similarity_score_repeats > 0.999:
                    self.similarity_score_repeats = 1.0
                try:
                    os.remove('file_for_similarity_{}.fa'.format(self.index))
                except OSError:
                    pass
            else:
                self.similarity_score_repeats = 1.0
        else:
            self.similarity_score_repeats = 1.0

        if self.list_spacers:
            if self.check_list(self.list_spacers):
                self._write_input_file_for_similarity(self.list_spacers)
                self.similarity_score_spacers = self._compute_similarity()
                try:
                    os.remove('file_for_similarity_{}.fa'.format(self.index))
                except OSError:
                    pass
            else:
                self.similarity_score_spacers = 1.0
        else:
            self.similarity_score_spacers = 1.0

    def output(self):
        return float(self.similarity_score_repeats), float(self.similarity_score_spacers)


class CrisprSimilarityNew(RepeatSpacersFeatures):
    def __init__(self, index, list_repeats, list_spacers):

        super().__init__(index, list_repeats, list_spacers)

        if self.list_spacers:
            if self.list_spacers[-1] == '':
                self.list_spacers = self.list_spacers[:-1]

        if set(self.list_spacers) == {"-"}:
            self.list_spacers = []

        self.similarity_score_repeats = None
        self.similarity_score_spacers = None
        self._compute_similarity_repeats_spacers()

    @staticmethod
    def check_list(input_list):
        for element in input_list:
            if len(set(element) - {"", "-"}) > 0:
                return True
        return False

    def _compute_similarity_repeats(self):
        seqs = []
        for index, repeat in enumerate(self.list_repeats):
            if repeat:
                seqs.append((str(index), repeat))
        graphs = list(seq2graph(seqs))
        x = kernel_matrix(graphs, r=3, d=4)
        n = x.shape[0]

        if n == 1:
            return 1.0
        else:
            return (np.sum(x) - n) / (n ** 2 - n)

    def _compute_similarity_spacers(self):
        seqs = []
        for index, spacer in enumerate(self.list_spacers):
            if spacer:
                seqs.append((str(index), spacer))

        graphs = list(seq2graph(seqs))
        x = kernel_matrix(graphs, r=3, d=4)

        n = x.shape[0]

        if n == 1:
            return 1.0
        else:
            return (np.sum(x) - n) / (n ** 2 - n)

    def _compute_similarity_repeats_spacers(self):
        if self.list_repeats:
            if self.check_list(self.list_repeats):
                self.similarity_score_repeats = self._compute_similarity_repeats()
                if self.similarity_score_repeats > 0.999:
                    self.similarity_score_repeats = 1.0
            else:
                self.similarity_score_repeats = 1.0
        else:
            self.similarity_score_repeats = 1.0

        if self.list_spacers:
            if self.check_list(self.list_spacers):
                self.similarity_score_spacers = self._compute_similarity_spacers()
            else:
                self.similarity_score_spacers = 1
        else:
            self.similarity_score_spacers = 1

    def output(self):
        return float(self.similarity_score_repeats), float(self.similarity_score_spacers)


class CrisprHmmer(RepeatSpacersFeatures):
    def __init__(self, index, list_repeats, list_spacers):

        super().__init__(index, list_repeats, list_spacers)
        self.crispr_seq = ''

        self.hmm_score = None

        self._create_crispr_seq()
        self._create_file()
        self._call_prodigal()
        self._run_hmm_search('../tools/hmm_search/models_tandem.hmm')
        self._get_best_score_from_hmm()
        self._clean_up()

    def _create_crispr_seq(self):
        for r, s in zip(self.list_repeats, self.list_spacers):
            self.crispr_seq += (r + s)
        self.crispr_seq += self.list_repeats[-1]

    def _create_file(self):
        with open('fasta_to_do_hmm_{}.fa'.format(self.index), 'w') as f:
            f.write('>seq_to_score\n')
            f.write(self.crispr_seq)

    def _call_prodigal(self):
        #cmd = 'tools/prodigal/prodigal'
        #cmd += ' -i fasta_to_do_hmm_{}.fa -p meta -a protein_{}.fa'.format(self.index, self.index)

        no_binary_cmd = f"prodigal -i fasta_to_do_hmm_{self.index}.fa -p meta -a protein_{self.index}.fa"
        process = subprocess.Popen(no_binary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        a, b = process.communicate()

    def _run_hmm_search(self, hmm_model):
        if os.stat('protein_{}.fa'.format(self.index)).st_size != 0:
            #cmd = 'tools/hmm_search/hmmsearch --tblout {} {} {}'.format('result_hmm_{}.out'.format(self.index),
            #                                                            hmm_model,
            #                                                            'protein_{}.fa'.format(self.index))

            no_binary_cmd = f"hmmsearch --tblout result_hmm_{self.index}.out {hmm_model} protein_{self.index}.fa"

            process = subprocess.Popen(no_binary_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            process.communicate()
        else:
            self.hmm_score = 0.0

    def _get_best_score_from_hmm(self):
        if self.hmm_score is None:
            with open('result_hmm_{}.out'.format(self.index), 'r') as f:
                list_lines = f.readlines()

            list_scores = []

            information_lines = list_lines[3:]
            for line in information_lines:
                list_info = line.split()
                if len(list_info) > 6:
                    score = float(list_info[5])
                    list_scores.append(score)
                else:
                    break

            self.hmm_score = max(list_scores) if list_scores else 0.0

    def _clean_up(self):
        try:
            os.remove('result_hmm_{}.out'.format(self.index))
        except OSError:
            pass
        finally:
            os.remove('protein_{}.fa'.format(self.index))
            os.remove('fasta_to_do_hmm_{}.fa'.format(self.index))

    def output(self):
        return float(self.hmm_score)


class FeatureExtractor:
    def __init__(self, index, crispr_candidate, features_to_extract=None):
        self.index = index
        self.crispr_candidate = crispr_candidate
        self.features_to_extract = features_to_extract
        self.list_repeats = crispr_candidate.list_repeats
        self.list_repeats_gaped = crispr_candidate.list_repeats_gaped
        self.list_spacers = crispr_candidate.list_spacers

        self.dict_features = {}

    def extract(self):
        gaped_consensus = self.crispr_candidate.consensus_gaped

        if 'number_repeats' in self.features_to_extract:
            number_repeats = CrisprNumberRepeats(self.index, self.list_repeats).output()
            self.dict_features['number_repeats'] = number_repeats

        if 'repeat_len' in self.features_to_extract:
            repeat_len = CrisprRepeatLen(self.index, self.list_repeats).output()
            self.dict_features['repeat_len'] = repeat_len

        if 'number_mismatches' in self.features_to_extract:
            number_mismatches = CrisprNumberMismatches(self.index, self.list_repeats_gaped, gaped_consensus).output()
            self.dict_features['number_mismatches'] = number_mismatches

        if 'avg_spacer_len' in self.features_to_extract:
            avg_spacer_len = CrisprAvgSpacerLength(self.index, self.list_spacers).output()
            self.dict_features['avg_spacer_len'] = avg_spacer_len

        if 'at_richness' in self.features_to_extract:
            at_richness = CrisprATRich(self.index, self.list_repeats, self.list_spacers).output()
            self.dict_features['at_richness'] = at_richness

        if 'spacer_evenness' in self.features_to_extract:
            spacer_eveness = CrisprSpacerEveness(self.index, self.list_spacers).output()
            self.dict_features['spacer_evenness'] = spacer_eveness

        if ('repeat_similarity' in self.features_to_extract) or ('spacer_similarity' in self.features_to_extract):
            repeat_similarity, spacer_similarity = CrisprSimilarityNew(self.index, self.list_repeats,
                                                                       self.list_spacers).output()
            self.dict_features['repeat_similarity'] = repeat_similarity
            self.dict_features['spacer_similarity'] = spacer_similarity

        if 'hmmr_score' in self.features_to_extract:
            hmmr_score = CrisprHmmer(self.index, self.list_repeats, self.list_spacers).output()
            self.dict_features['hmmr_score'] = hmmr_score

        feature_vector = [self.dict_features[feature] for feature in self.features_to_extract]

        return np.asarray(feature_vector).reshape(1, -1)



#                 HELPERS FOR SCORING
##############################################################
##############################################################


def get_full_vector(list_vectors):
    """
    8: (2, 4, 5, 6, 7, 8, 9, 11),
    9: (1, 2, 4, 5, 7, 8, 9, 10, 12),
    10: (0, 2, 3, 4, 5, 6, 7, 10, 11, 12),
    """
    fist_model_vector = list_vectors[0][0]
    second_model_vector = list_vectors[1][0]
    third_model_vector = list_vectors[2][0]

    f_0 = third_model_vector[0]
    f_1 = second_model_vector[0]
    f_2 = fist_model_vector[0]
    f_3 = third_model_vector[2]
    f_4 = third_model_vector[3]
    f_5 = third_model_vector[4]
    f_6 = third_model_vector[5]
    f_7 = third_model_vector[6]
    f_8 = fist_model_vector[5]
    f_9 = second_model_vector[6]
    f_10 = second_model_vector[7]
    f_11 = fist_model_vector[7]
    f_12 = second_model_vector[8]

    return np.asarray([f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11, f_12]).reshape(1, -1)