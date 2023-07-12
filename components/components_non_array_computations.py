import os
import subprocess
from os import listdir
from os.path import isfile, join
import shutil
import csv

from components.components_detection_refinement import CrisprCandidate

#              Strand computation
#######################################################
#######################################################


def rev_compliment(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', "K": "K"}
    try:
        compliment_seq = "".join([complement[nt] for nt in seq])
    except KeyError:
        compliment_seq = ""
        for char in seq:
            if char in complement:
                compliment_seq += complement[char]
            else:
                compliment_seq += char
    return compliment_seq[::-1]


def to_rna(seq):
    return seq.replace("T", "U")


def get_orientation(sequence, absolute_directory_path):
    rev_comp = rev_compliment(sequence)

    orig_rna, rev_comp_rna = to_rna(sequence), to_rna(rev_comp)

    with open("forward.txt", "w") as f:
        f.write(orig_rna)

    with open("reversed.txt", "w") as f:
        f.write(rev_comp_rna)

    cmd = f"{absolute_directory_path}/tools/strand_prediction/EDeN -a TEST -i forward.txt -M 1 -r 3 -d 3 -f SEQUENCE -g" \
          " DIRECTED -m tools/strand_prediction/DR_Repeat_model"

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.communicate()
    f = open("prediction", "r")
    val_in_plus = float((f.readline().split())[1])
    f.close()
    cmd = f"{absolute_directory_path}/tools/strand_prediction/EDeN -a TEST -i reversed.txt -M 1 -r 3 -d 3 -f SEQUENCE -g" \
          " DIRECTED -m tools/strand_prediction/DR_Repeat_model"

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    process.communicate()
    f = open("prediction", "r")
    val_in_minus = float((f.readline().split())[1])
    f.close()

    os.remove("forward.txt")
    os.remove("reversed.txt")
    os.remove("prediction")
    val_in_minus *= -1

    val = (val_in_plus + val_in_minus) / 2
    if val > 0:
        return 1
    else:
        return 0


class StrandComputationSingleCrispr:
    def __init__(self, crispr_array, absolute_directory_path):
        self.crispr_array = crispr_array
        self.absolute_directory_path = absolute_directory_path
        self.strand = None

        self._compute_strand()

    def _compute_strand(self):
        consensus = self.crispr_array.consensus
        self.strand = get_orientation(consensus, self.absolute_directory_path)
        self.strand = "Forward" if self.strand else "Reversed"

    def output(self):
        return self.strand


class StrandComputation:
    def __init__(self, list_of_crisprs, absolute_directory_path):
        self.list_of_crisprs = list_of_crisprs
        self.absolute_directory_path = absolute_directory_path
        self.dict_strands = {}

        self._compute_all_strands()

    def _compute_all_strands(self):
        for index, crispr in enumerate(self.list_of_crisprs):
            consensus = crispr.consensus
            strand = get_orientation(consensus, self.absolute_directory_path)
            strand = "Forward" if strand else "Reversed"
            self.dict_strands[index] = strand

    def output(self):
        return self.dict_strands


class StrandComputationNew:
    def __init__(self, list_of_crisprs, absolute_directory_path):
        self.list_of_crisprs = list_of_crisprs
        self.absolute_directory_path = absolute_directory_path
        self.dict_strands = {}

        self._compute_all_strands()

    def _compute_all_strands(self):
        if self.list_of_crisprs:
            with open("CRISPR_arrays_for_strand.fa", "w") as f:
                for array_index, array in enumerate(self.list_of_crisprs, 1):
                    consensus = array.consensus
                    consensus_fidex_alphabet = self.remove_non_canonical_char_from_string(consensus)
                    f.write(f">CRISPR_{array_index}_consensus\n{consensus_fidex_alphabet}\n")

                if len(self.list_of_crisprs) == 2:
                    consensus = self.list_of_crisprs[-1].consensus
                    consensus_fidex_alphabet = self.remove_non_canonical_char_from_string(consensus)
                    f.write(f">CRISPR_3_consensus\n{consensus_fidex_alphabet}\n")

            try:
                os.mkdir("ResultsStrand")
            except Exception:
                pass

            cmd = f"python {self.absolute_directory_path}/tools/strand_prediction/CRISPRstrand/CRISPRstrand.py -r -i CRISPR_arrays_for_strand.fa --model_path {self.absolute_directory_path}/tools/strand_prediction/CRISPRstrand/Models/model_r.h5 --output_folder ResultsStrand"
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            a, b = process.communicate()
            #print(a, b)
            #os.remove("CRISPR_arrays_for_strand.fa")

            with open("ResultsStrand/CRISPRstrand_Summary.tsv", "r") as f:
                lines = f.readlines()
                strands = [l.split()[3] for l in lines[1:]]
                strands = ["Reversed" if s == "Reverse" else "Forward" for s in strands]

            for index, crispr in enumerate(self.list_of_crisprs):
                self.dict_strands[index] = strands[index]

            try:
                shutil.rmtree("ResultsStrand")
            except Exception:
                pass

    @staticmethod
    def remove_non_canonical_char_from_string(consensus_string):
        canonical_chars = ["A", "T", "G", "C"]
        new_consensus_string = ""
        for char in consensus_string:
            if char in canonical_chars:
                new_consensus_string += char
            else:
                new_consensus_string += "N"
        return new_consensus_string

    def output(self):
        return self.dict_strands



#      IS element computation
################################################
################################################


class FastaMatch:
    def __init__(self, first_id, second_id, similarity, coverage):
        self.first_id = int(first_id)
        self.second_id = int(second_id)
        self.similarity = similarity
        self.coverage = coverage


class FastaSimilarity():
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


#          Components cas genes
######################################################
######################################################

def parse_csv_protein_file(file_name):
    dict_file_csv = {}
    with open(file_name) as f:
        lines = f.readlines()[1:]
    for line in lines:
        line_info = line.split(",")
        start = int(line_info[1])
        end = int(line_info[2])
        annotation = line_info[5].strip()
        key = (start, end)
        if "cas" in annotation:
            dict_file_csv[key] = annotation

    return dict_file_csv


def cas_identifier_result_folder_parser(folder_path):
    dict_cas_proteins = {}
    onlyfiles = [f for f in listdir(folder_path) if isfile(join(folder_path, f))]
    protein_files = [f for f in onlyfiles if "annotated_proteins" in f]
    for file_name in protein_files:
        full_path = join(folder_path, file_name)
        dict_file = parse_csv_protein_file(full_path)
        dict_cas_proteins = {**dict_cas_proteins, **dict_file}
    return dict_cas_proteins


def cas_identifier_cassete_csv_parser(file_path):
    cassette_dict = {}
    if isfile(file_path):
        with open(file_path) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            next(csv_reader)  # skip header row
            for row in csv_reader:
                cassette_id = int(row[6])
                start = int(row[1])
                end = int(row[2])
                if cassette_id not in cassette_dict:
                    cassette_dict[cassette_id] = (start, end)
                else:
                    current_start, current_end = cassette_dict[cassette_id]
                    cassette_dict[cassette_id] = (min(current_start, start), max(current_end, end))
    return cassette_dict


def cas_identifier_read_predicted_labels(file_path):
    dict_predicted_labels = {}
    if isfile(file_path):
        with open(file_path) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            next(csv_reader)  # skip header row
            for row in csv_reader:
                cassette_id = int(row[1])
                predicted_label = row[4]
                dict_predicted_labels[cassette_id] = predicted_label
    return dict_predicted_labels


def cas_identifier_combine_dicts(cassette_dict, predicted_labels):
    combined_dict = {}
    for cassette_id in cassette_dict:
        if cassette_id in predicted_labels:
            start, end = cassette_dict[cassette_id]
            label = predicted_labels[cassette_id]
            combined_dict[cassette_id] = (start, end, label)
    return combined_dict


def run_cas_identifier(file_name, absolute_directory_path):
    try:
        #cmd1 = f"mkdir {absolute_directory_path}/output_cas"
        #cmd2 = f"mkdir {absolute_directory_path}/output_cas/cassette"

        cmd1 = f"mkdir output_cas"
        cmd2 = f"mkdir output_cas/cassette"

        process = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        process = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

    except Exception:
        pass

    try:
        #cmd = f"mkdir {absolute_directory_path}/output"
        cmd = f"mkdir output"
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()
    except Exception:
        pass

    command = f"python {absolute_directory_path}/tools/CRISPRcasIdentifier/CRISPRcasIdentifier/CRISPRcasIdentifier.py -f {file_name} -ho output_cas/hmmsearch -st dna -co output_cas/cassette"
    #print(command)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    a, b = process.communicate()
    #print(a, b)

    directory_path = os.path.dirname(file_name)
    log_file_name_log = os.path.splitext(os.path.basename(file_name))[0] + "_prodigal.log"
    log_file_path_log = os.path.join(directory_path, log_file_name_log)

    log_file_name_fa = os.path.splitext(os.path.basename(file_name))[0] + "_proteins.fa"
    log_file_path_fa = os.path.join(directory_path, log_file_name_fa)

    #print(log_file_path_log)
    #print(log_file_path_fa)

    #folder = "/".join(file_name.split("/")[:-1])
    #file_base = file_name.split("/")[-1].split(".")[0]
    #log_prodigal_file = file_base + "_prodigal.log"
    #prodigal_prots = file_base + "_proteins.fa"
    #full_path_log = join(folder, log_prodigal_file)
    #full_path_prots = join(folder, prodigal_prots)

    try:
        os.remove(log_file_path_log)
    except Exception:
        pass

    try:
        os.remove(log_file_path_fa)
    except Exception:
        pass


def complete_info_with_cas_identifier(file_name, absolute_directory_path):
    run_cas_identifier(file_name, absolute_directory_path)
    #cas_path = f"{absolute_directory_path}/output_cas/cassette"
    #cas_path_short = f"{absolute_directory_path}/output_cas/"

    cas_path = f"output_cas/cassette"
    cas_path_short = f"output_cas/"

    dict_cas = cas_identifier_result_folder_parser(cas_path)
    dict_casette_classes = cas_identifier_read_predicted_labels("output/predictions.csv")
    dict_cassete_intervals = cas_identifier_cassete_csv_parser(f"{cas_path}/HMM2019_cassettes.csv")
    dict_casete_start_end_class = cas_identifier_combine_dicts(dict_cassete_intervals, dict_casette_classes)
    try:
        shutil.rmtree(cas_path)
    except Exception:
        pass

    try:
        shutil.rmtree(cas_path_short)
    except Exception:
        pass



    return dict_cas, dict_casete_start_end_class


#        Leader seq search
##########################################################
##########################################################


def rev_compliment_seq(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', ' ': ' ', "-": "-", ".": "."}
    try:
        compliment_seq = "".join([complement[nt] for nt in seq])
    except KeyError:
        compliment_seq = ""
        for nt in seq:
            if nt in complement:
                compliment_seq += complement[nt]
            else:
                compliment_seq += nt
    return compliment_seq[::-1]


class FullLeaderSeqSearch:
    def __init__(self, list_crispr_candidates, list_corresponding_strands, full_dna):
        self.list_crispr_candidates = list_crispr_candidates
        self.list_strands = list_corresponding_strands
        self.full_dna = full_dna

        self.dict_leader_data = {}
        self.dict_downstream_data = {}

        self._compute_all_leaders()

    def _compute_all_leaders(self):
        for index, strand, crispr_candidate in zip(range(len(self.list_crispr_candidates)),
                                                   self.list_strands,
                                                   self.list_crispr_candidates):
            leader, downstream = LeaderSeqSearch(crispr_candidate, strand, self.full_dna).output()
            self.dict_leader_data[index] = leader
            self.dict_downstream_data[index] = downstream

    def output(self):
        return self.dict_leader_data, self.dict_downstream_data


class LeaderSeqSearch:
    def __init__(self, crispr_candidate, strand, full_dna):
        self.crispr_candidate = crispr_candidate
        self.strand = strand
        self.full_dna = full_dna

        self.leader = None
        self.downstream_region = None

        self._compute_leader_seq()

    def _compute_leader_seq(self):
        list_repeat_indexes = self.crispr_candidate.list_repeat_starts
        list_repeats = self.crispr_candidate.list_repeats
        first_index = list_repeat_indexes[0]

        array_end = list_repeat_indexes[-1]
        last_repeat_len = len(list_repeats[-1])

        if self.strand == "Forward":
            start = first_index - 201
            end = first_index - 1
            self.leader = self.full_dna[start:end+1]
            self.downstream_region = self.full_dna[array_end+last_repeat_len:array_end+last_repeat_len+200]
        else:
            start = first_index - 201
            end = first_index - 1
            self.downstream_region = self.full_dna[start:end]
            self.leader = self.full_dna[array_end + last_repeat_len:array_end + last_repeat_len + 200]

            self.downstream_region = rev_compliment_seq(self.downstream_region)
            self.leader = rev_compliment_seq(self.leader)

    def output(self):
        return self.leader, self.downstream_region

#           RevCom computations
#################################################
#################################################


class RevComComputation:
    def __init__(self, crispr_candidate):
        self.forward_crispr_candidate = crispr_candidate
        self.rev_com_candidate = None
        self._compute_rev_com_candidate()

    def _compute_rev_com_candidate(self):
        list_repeats_f = self.forward_crispr_candidate.list_repeats
        list_repeats_gaped_f = self.forward_crispr_candidate.list_repeats_gaped
        list_repeat_starts = self.forward_crispr_candidate.list_repeat_starts
        list_spacers = self.forward_crispr_candidate.list_spacers

        list_repeats_rev_com = [rev_compliment_seq(repeat) for repeat in list_repeats_f][::-1]

        list_repeats_gaped_rev_com = [rev_compliment_seq(repeat_gaped)
                                      for repeat_gaped in list_repeats_gaped_f][::-1]

        list_repeat_starts_rev_com = list_repeat_starts[::-1]

        list_spacers_rev_com = [spacer for spacer in list_spacers][::-1]

        self.rev_com_candidate = CrisprCandidate(list_repeats=list_repeats_rev_com,
                                                 list_repeats_gaped=list_repeats_gaped_rev_com,
                                                 list_repeat_starts=list_repeat_starts_rev_com,
                                                 list_spacers=list_spacers_rev_com)

    def output(self):
        return self.rev_com_candidate
