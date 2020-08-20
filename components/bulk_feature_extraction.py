import numpy as np
from collections import namedtuple
import math
import os
import re
from time import time
import subprocess

import sys
import joblib

sys.path.insert(0, 'components/')
sys.path.insert(0, 'tools/')


class BulkFeatureExtractorOrf:
    def __init__(self, dict_crispr_candidates):
        self.dict_crispr_candidates = dict_crispr_candidates

        #print("Orf")
        start_time = time()

        self._create_repeat_intervals_for_ofr()
        self._extract_orf_scores()
        self._clean_up()
        end_time = time()
        time_taken = end_time - start_time
        #print(f"Orf took {time_taken}")

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

        cmd = 'tools/prodigal/prodigal'
        cmd += ' -i file_for_prodigal.fa -o prodigal_result.txt -c -p meta'
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
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
            pass
            #os.remove("file_for_prodigal.fa")
        except OSError:
            pass

        try:
            pass
            #os.remove("prodigal_result.txt")
        except OSError:
            pass

    def output_ofr_scores(self):
        return self.dict_final_orf_result


class BulkFeatureExtractorORFHMMR:
    def __init__(self, dict_crispr_candidates):
        self.dict_crispr_candidates = dict_crispr_candidates

        start_time = time()

        self._create_repeat_intervals_for_ofr()
        self._extract_orf_scores_and_proteins()
        self._run_hmm_search()
        self._extract_best_score_from_hmm()
        self._clean_up()
        end_time = time()
        time_taken = end_time - start_time
        # print(f"Orf took {time_taken}")

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

                    #f.write(f'>{key}_{index}\n')
                    start, end = key
                    f.write(f">{start}_{end}_{index}\n")
                    f.write(crispr_seq)
                    f.write("\n")

        cmd = 'tools/prodigal/prodigal'
        cmd += ' -i file_for_prodigal.fa -o prodigal_result.txt -c -p meta'
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        cmd = 'tools/prodigal/prodigal'
        cmd += ' -i file_for_prodigal.fa -p meta -a protein_results.fa'
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
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
        cmd = 'tools/hmm_search/hmmsearch --tblout result_hmm.out tools/hmm_search/models_tandem.hmm protein_results.fa'
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
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

        #print("Blast")
        start_time = time()

        self.dict_blast_scores_1 = {}
        self.dict_blast_scores_2 = {}
        self._extract_all_blast_scores()
        self._clean_up()

        end_time = time()
        time_taken = end_time - start_time
        #print(f"Blast took {time_taken}")

    def _extract_all_blast_scores(self):
        list_consensus_rep = []
        for key, list_crispr_candidates in self.dict_crispr_candidates.items():
            for index, crispr_candidate in enumerate(list_crispr_candidates):
                list_consensus_rep.append(crispr_candidate.consensus)

        start_time = time()

        with open('file_with_all_consensus.fa', 'w') as f:
            for index, consensus_seq in enumerate(list_consensus_rep, 1):
                f.write(f">{index}\n")
                f.write(f"{consensus_seq}\n")

        db_file = 'Verified_repeats_dataset1.fa'

        cmd = 'tools/blasting/blastn -query file_with_all_consensus.fa'
        cmd += ' -db tools/blasting/'
        cmd += db_file
        cmd += ' -word_size 6'
        cmd += ' -outfmt 6 -out output_fasta_bulk_extraction1'

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        db_file = 'Verified_repeats_dataset2.fa'

        cmd = 'tools/blasting/blastn -query file_with_all_consensus.fa'
        cmd += ' -db tools/blasting/'
        cmd += db_file
        cmd += ' -word_size 6'
        cmd += ' -outfmt 6 -out output_fasta_bulk_extraction2'

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        end_time = time()
        only_blast = end_time - start_time
        #print(f"Only blast tool {only_blast}")

        with open("output_fasta_bulk_extraction1") as f:
            lines = f.readlines()
            #max_num = int(lines[-1].split()[0])
            max_num = len(list_consensus_rep)
            first_lines_scores = {}

            for line in lines:
                input_index = int(line.split()[0])
                if input_index not in first_lines_scores:
                    first_lines_scores[input_index] = float(line.strip().split()[-1])

            scores1 = []

            for possible_input_input_index in range(1, max_num+1):
                #print(possible_input_input_index)
                if possible_input_input_index in first_lines_scores:
                    bit_score = first_lines_scores[possible_input_input_index]
                else:
                    bit_score = 0.0

                scores1.append(bit_score)

        with open("output_fasta_bulk_extraction2") as f:
            lines = f.readlines()
            #max_num = int(lines[-1].split()[0])
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
        #print("MFE")
        start_time = time()

        self.dict_crispr_candidates = dict_crispr_candidates
        self.dict_final_mfe_result = {}

        self._extract_mfe_scores()
        self._clean_up()

        end_time = time()
        time_taken = end_time - start_time
        #print(f"MFE took {time_taken}")

    def _extract_mfe_scores(self):

        with open("file_for_mfe.fa", "w") as f:
            for key, list_crispr_candidates in self.dict_crispr_candidates.items():
                for index, crispr_candidate in enumerate(list_crispr_candidates):
                    consensus = crispr_candidate.consensus

                    f.write(f'>{key}_{index}\n')
                    f.write(consensus)
                    f.write("\n")

        cmd = "cat file_for_mfe.fa | tools/rna_fold/RNAfold --noLP --noPS > rna_fold_output.txt"

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
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
