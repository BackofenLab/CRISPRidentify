import os
import collections
from collections import namedtuple
import math
import os
import subprocess
import re


class CRTParser:
    def __init__(self, crt_file_name):
        self.crt_file_name = crt_file_name
        self.lines_in_file = []
        self.list_crisprs = []

        self.get_the_file()
        self.parse_the_lines()

    def get_the_file(self):
        with open(self.crt_file_name) as f:
            self.lines_in_file = f.readlines()

    def parse_the_lines(self):
        flag_in = False
        index = 1
        list_lines_crispr = []
        for line in self.lines_in_file:
            if '--------' in line:
                if flag_in:
                    crispr = CRTCrispr(list_lines=list_lines_crispr)
                    self.list_crisprs.append(crispr)
                    index += 1
                    list_lines_crispr = []
                flag_in = not flag_in
            elif flag_in:
                list_lines_crispr.append(line)


class CRTCrispr:
    def __init__(self, list_lines):
        self.list_lines = list_lines
        self.list_repeats = []
        self.list_spacers = []
        self.consensus_rep = ''
        self.len_consensus = None
        self.cr_start = None
        self.cr_end = None

        self._get_repeats_spacers()
        self._get_consensus()
        self._get_location()

    def _get_repeats_spacers(self):
        for line in self.list_lines:
            rep = line.split('\t')[2]
            sp = line.split('\t')[3]
            self.list_repeats.append(rep)
            if len(sp.strip()) > 0:
                self.list_spacers.append(sp)

    def _get_consensus(self):
        for char_ind in range(len(self.list_repeats[0])):
            list_char_with_index = [repeat[char_ind] for repeat in self.list_repeats]
            most_common_char = max(set(list_char_with_index), key=list_char_with_index.count)
            self.consensus_rep += most_common_char
        self.len_consensus = len(self.consensus_rep)

    def _get_location(self):
        self.cr_start = int(self.list_lines[0].split('\t')[0])
        self.cr_end = int(self.list_lines[-1].split('\t')[0]) + len(self.list_repeats[-1]) - 1


class FastaMatch:
    def __init__(self, first_id, second_id, similarity, coverage):
        self.first_id = int(first_id)
        self.second_id = int(second_id)
        self.similarity = similarity
        self.coverage = coverage


class FastaSimilarity(object):
    def __init__(self, list_sequences, similarity, coverage):
        self.list_sequences = list_sequences
        self.similarity = similarity
        self.coverage = coverage

        self.fasta_info = None
        self.pairs = []

        self._run_fasta()
        self._parse_fasta()

        self._filter_by_similarity()
        self._filter_by_coverage()

    def _run_fasta(self):
        with open("sequences.fa", "w") as f:
            for index, repeat in enumerate(self.list_sequences):
                f.write(">{}\n".format(index))
                f.write("{}\n".format(repeat.strip()))
        cmd = "tools/fasta/fasta36 sequences.fa sequences.fa -m 8 > fasta_similarity.fastab"

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()
        os.remove("sequences.fa")

    def _parse_fasta(self):
        with open("fasta_similarity.fastab", "r") as f:
            lines = f.readlines()

        self.fasta_info = [line.strip().split("\t") for line in lines]
        self.fasta_info = [fasta_line_match for fasta_line_match in self.fasta_info
                           if int(fasta_line_match[0]) < int(fasta_line_match[1])]

        os.remove("fasta_similarity.fastab")

    def _filter_by_similarity(self):
        self.fasta_info = [fasta_line_match for fasta_line_match in self.fasta_info
                           if float(fasta_line_match[2]) / 100 >= self.similarity - 0.00005]

    def _filter_by_coverage(self):
        filtered_by_coverage = []
        for fasta_line_match in self.fasta_info:
            interval_first_element = abs(int(fasta_line_match[6]) - int(fasta_line_match[7])) + 1
            interval_second_element = abs(int(fasta_line_match[8]) - int(fasta_line_match[9])) + 1

            coverage_first = float(interval_first_element) / len(self.list_sequences[int(fasta_line_match[0])])
            coverage_second = float(interval_second_element) / len(self.list_sequences[int(fasta_line_match[1])])
            min_coverage = min(coverage_first, coverage_second)
            if min_coverage >= self.coverage - 0.00005:
                min_coverage = int(min_coverage * 100)
                fasta_line_match.append(min_coverage)
                filtered_by_coverage.append(fasta_line_match)

        self.fasta_info = filtered_by_coverage

    def output(self):
        list_fasta_info_objects = []
        for line in self.fasta_info:
            fasta_line_match_important = [line[i] for i in (0, 1, 2, -1)]
            fasta_match = FastaMatch(*fasta_line_match_important)
            list_fasta_info_objects.append(fasta_match)
        return list_fasta_info_objects


class HMMMatch:
    def __init__(self, target, query, e_value, score):
        self.target = target
        self.query = query
        self.e_value = e_value
        self.score = score

    def __repr__(self):
        return "{} {} {} {}".format(self.target, self.query, self.e_value, self.score)


class HMMMatchProteinCoordinates:
    def __init__(self, target, query,  e_value, score, protein_start, protein_end, protein_strand):
        self.target = target
        self.query = query
        self.e_value = e_value
        self.score = score
        self.protein_start = protein_start
        self.protein_end = protein_end
        self.protein_strand = protein_strand

    def __repr__(self):
        return "{} {} {} {} {} {} {}".format(self.target, self.query, self.e_value,
                                             self.score, self.protein_start,
                                             self.protein_end, self.protein_strand)


class HMMResultParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.list_matches = []
        self._input_parsing()

    def _input_parsing(self):
        with open(self.file_path, "r") as f:
            information_lines = f.readlines()

        information_lines = information_lines[3:]
        for line in information_lines:
            list_info = line.split()
            if list_info[0] != "#":
                score = float(list_info[5])
                e_value = float(list_info[4])
                target = list_info[0]
                query = list_info[2]

                hmm_match = HMMMatch(target, query, e_value, score)
                self.list_matches.append(hmm_match)

    def output(self):
        return self.list_matches


class CrisprPair:
    def __init__(self, first_interval, second_interval, first_consensus,
                 second_consensus, similarity, coverage, first_strand="R", second_strand="R"):

        self.first_interval = first_interval
        self.second_interval = second_interval
        self.first_consensus = first_consensus
        self.second_consensus = second_consensus
        self.similarity = similarity
        self.coverage = coverage
        self.first_strand = first_strand
        self.second_strand = second_strand

    def __repr__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(str(self.first_interval), self.first_consensus,
                                                       str(self.second_interval), self.second_consensus,
                                                       self.similarity, self.coverage, self.first_strand,
                                                       self.second_strand)

    def coordinates_between(self):
        return self.first_interval[1], self.second_interval[0]


class ISElementSearch:
    def __init__(self, dna_seq, hmm_model):
        self.dna_seq = dna_seq
        self.hmm_model = hmm_model

        self.dict_protein_description = {}
        self.hmm_matches = []
        self.hmm_protein_match = []

        self._protein_extraction()
        self._hmm_extraction()
        self._hmm_match_protein_creation()

    def _protein_extraction(self):
        with open('dna.fa', 'w') as f:
            f.write('>seq_to_score\n')
            f.write(self.dna_seq)

        cmd = 'tools/prodigal/prodigal'
        cmd += ' -i dna.fa -p meta -a protein.fa'
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()
        os.remove('dna.fa')

        with open("protein.fa", "r") as f:
            lines = f.readlines()
        header_lines = [line[1:] for line in lines if ">" in line]
        keys = [line.split()[0] for line in header_lines]
        coordinates = [(int(line.split()[2]), int(line.split()[4])) for line in header_lines]
        strands = [int(line.split()[6]) for line in header_lines]
        strands = ["Forward" if x == 1 else "Reverse" for x in strands]
        values = [(*coordinate, strand) for coordinate, strand in zip(coordinates, strands)]

        self.dict_protein_description = {key: value for key, value in zip(keys, values)}

    def _hmm_extraction(self):
        if os.stat('protein.fa').st_size != 0:
            cmd = 'tools/hmm_search/hmmsearch --tblout {} {} {}'.format('result_hmm.out', self.hmm_model,
                                                                        'protein.fa')
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            process.communicate()

            hmm_result_parser = HMMResultParser("result_hmm.out")
            self.hmm_matches = hmm_result_parser.output()
            os.remove('result_hmm.out')

        os.remove('protein.fa')

    def _hmm_match_protein_creation(self):
        for hmm_match in self.hmm_matches:
            hmm_target = hmm_match.target
            hmm_query = hmm_match.query
            hmm_e_value = hmm_match.e_value
            hmm_score = hmm_match.score
            protein_start, protein_end, protein_strand = self.dict_protein_description[hmm_target]
            self.hmm_protein_match.append(HMMMatchProteinCoordinates(hmm_target, hmm_query, hmm_e_value,
                                                                     hmm_score, protein_start, protein_end,
                                                                     protein_strand))

    def output(self):
        return self.hmm_protein_match


class FullISElementSearch:
    def __init__(self, full_dna, list_of_crisprs, hmm_model, min_similarity, min_coverage):
        self.full_dna = full_dna
        self.list_of_crisprs = list_of_crisprs
        self.hmm_model = hmm_model
        self.min_similarity = min_similarity
        self.min_coverage = min_coverage

        self.crispr_pairs = []
        self.is_results = {}

        self._compute_all_is_elements()

    def _compute_all_is_elements(self):
        self._get_consensus_repeats()
        self._get_crispr_pairs()
        self._perform_hmm_search()

    def _get_consensus_repeats(self):
        self.consensus_repeats = [crispr.consensus for crispr in self.list_of_crisprs]

    def _get_crispr_pairs(self):
        fasta_similarity = FastaSimilarity(self.consensus_repeats, self.min_similarity, self.min_coverage)
        fasta_results = fasta_similarity.output()
        for fasta_match in fasta_results:
            first_id, second_id = fasta_match.first_id, fasta_match.second_id
            if second_id - first_id == 1:
                cr_pair = (first_id, second_id)
                self.crispr_pairs.append(cr_pair)

    def _perform_hmm_search(self):
        for cr_pair in self.crispr_pairs:
            first_index, second_index = cr_pair[0], cr_pair[1]
            end_first = self.list_of_crisprs[first_index].list_repeat_starts[-1] + \
                        len(self.list_of_crisprs[first_index].list_repeats[-1])

            begin_second = self.list_of_crisprs[second_index].list_repeat_starts[0]

            dna_interval = self.full_dna[end_first:begin_second]
            is_element_search = ISElementSearch(dna_interval, self.hmm_model)
            is_search_result = is_element_search.output()
            sorted_is_search_result = sorted(is_search_result, key=lambda x: x.e_value)
            if sorted_is_search_result:
                best_hit = sorted_is_search_result[0]
                best_interval_start = best_hit.protein_start
                best_interval_end = best_hit.protein_end
                best_interval_strand = best_hit.protein_strand
                best_interval_target = best_hit.target
                best_interval_query = best_hit.query
                self.is_results[cr_pair[0]] = (end_first+best_interval_start, end_first+best_interval_end,
                                               best_interval_strand, best_interval_target, best_interval_query)

    def output(self):
        return self.is_results




