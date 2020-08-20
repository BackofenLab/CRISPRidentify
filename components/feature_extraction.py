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


from eden_rna import load
from eden_rna import fold as seq2graph
from eden_graph import kernel_matrix


class SaveLoad:
    def save(self, model_name):
        joblib.dump(self, model_name, compress=1)

    def load(self, model_name):
        self.__dict__.update(joblib.load(model_name).__dict__)


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






#class EdenEstimatorSaveLoad(EdenEstimator, SaveLoad):
#    pass

"""
class EdenClassifier:
    def __init__(self, load_option=None, hyper_parameters=None):
        self.hyper_parameters_eden = hyper_parameters
        self._load_option = load_option
        
        self._init_est()
        
    def _init_est(self):
        if self._load_option is not None:
            self._load_model()
        else:
            if not self.hyper_parameters_eden:
                self.est = EdenEstimatorSaveLoad(r=3, d=4, normalization=True, inner_normalization=True)
            else:
                self.est = EdenEstimatorSaveLoad(**self.hyper_parameters_eden)
                
    def _load_model(self):
        self.est = EdenEstimatorSaveLoad(r=3, d=4, normalization=True, inner_normalization=True)
        self.est.load(self._load_option)  
    
    def _make_data(self, pos_seqs, neg_seqs):
        pos_graphs = list(seq2graph(pos_seqs))
        neg_graphs = list(seq2graph(neg_seqs))

        graphs = pos_graphs + neg_graphs
        targets = np.array([1]*len(pos_graphs) + [0]*len(neg_graphs))
        return graphs, targets
        
    def train_eden_classifier(self, train_eden_pos, train_eden_neg):        
        graphs_train_eden, targets_train_eden = self._make_data(train_eden_pos, train_eden_neg)    
        self.est.fit(graphs_train_eden, targets_train_eden)            
    
    def test_eden_classifier(self, test_eden_pos, test_eden_neg):
        graphs_test, targets_test = self._make_data(test_eden_pos, test_eden_neg)
        preds = self.est.predict(graphs_test)        
        return 1 - np.count_nonzero(targets_test-preds)/float(len(preds))
    
    def score_input(self, seq_to_score):
        graphs = list(seq2graph(seq_to_score))
        scores = self.est.decision_function(graphs)
        return scores
    
    def save_model(self, model_name):
        self.est.save(model_name)

"""


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
        start_time = time()

        super().__init__(index, list_spacers)

        self.list_right_boundaries = list_right_boundaries
        
        self.spacer_eveness = None
        self._compute_spacer_eveness()

        end_time = time()
        time_taken = end_time - start_time
        #print(f"Spacer evenness took {time_taken}")
        
    def _compute_spacer_eveness(self):
        list_bin_values = [0 for _ in range((len(self.list_right_boundaries) + 1))]        
        for spacer in self.list_spacers:
            for boundary_index, right_boundary in enumerate(self.list_right_boundaries):
                if len(spacer) <= right_boundary:
                    list_bin_values[boundary_index] += 1
                    break
            else:
                list_bin_values[-1] += 1        

        list_bin_normalized_values = [x/float(sum(list_bin_values)) for x in list_bin_values]
        summ = sum([(-1) * x * math.log(x + 0.00000001) for x in list_bin_normalized_values])
        self.spacer_evenness = round(summ/(math.log(len(list_bin_values))), 6)
        
    def output(self):
        return float(self.spacer_evenness)


class CrisprSimilarity(RepeatSpacersFeatures):
    def __init__(self, index, list_repeats, list_spacers):
        start_time = time()

        super().__init__(index, list_repeats, list_spacers)

        if self.list_spacers:
            if self.list_spacers[-1] == '':
                self.list_spacers = self.list_spacers[:-1]

        if set(self.list_spacers) == {"-"}:
            self.list_spacers = []
        
        self.similarity_score_repeats = None
        self.similarity_score_spacers = None
        self._compute_similarity_repeats_spacers()

        end_time = time()
        time_taken = end_time - start_time
        #print(f"Similarity repeat_spacers took {time_taken}")

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
        # Compute the avg similarities between  sequences.
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
                self.similarity_score_spacers = 0
        else:
            self.similarity_score_spacers = 0
        
    def output(self):
        return float(self.similarity_score_repeats), float(self.similarity_score_spacers)


class CrisprSimilarityNew(RepeatSpacersFeatures):
    def __init__(self, index, list_repeats, list_spacers):
        start_time = time()

        super().__init__(index, list_repeats, list_spacers)

        if self.list_spacers:
            if self.list_spacers[-1] == '':
                self.list_spacers = self.list_spacers[:-1]

        if set(self.list_spacers) == {"-"}:
            self.list_spacers = []

        self.similarity_score_repeats = None
        self.similarity_score_spacers = None
        self._compute_similarity_repeats_spacers()

        end_time = time()
        time_taken = end_time - start_time
        #print(f"Similarity repeat_spacers took {time_taken}")

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

        st_time = time()
        graphs = list(seq2graph(seqs))
        x = kernel_matrix(graphs, r=3, d=4)
        end_time = time()
        #print("repeats", end_time - st_time)
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

        st_time = time()
        graphs = list(seq2graph(seqs))
        x = kernel_matrix(graphs, r=3, d=4)
        end_time = time()
        #print("spacers", end_time - st_time)

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
                self.similarity_score_spacers = 0
        else:
            self.similarity_score_spacers = 0

    def output(self):
        return float(self.similarity_score_repeats), float(self.similarity_score_spacers)


class CrisprORFScore(RepeatSpacersFeatures):
    def __init__(self, index, list_repeats, list_spacers):
        start_time = time()

        super().__init__(index, list_repeats, list_spacers)
        self.crispr_seq = ''
        self.dict_cds = {}
        self.list_repeat_intervals = []
        self.orf_best_score = 0.0
        
        self._create_repeat_intervals()
        self._create_crispr_seq()
        self._create_file()
        self._call_prodigal()
        self._read_from_orf_file()
        self._find_best_score()
        self._clean_up()

        end_time = time()
        time_taken = end_time - start_time
        #print(f"ORF score took {time_taken}")
        
    def _create_repeat_intervals(self):
        list_length_repeats = [len(repeat) for repeat in self.list_repeats]
        list_length_spacers = [len(spacer) for spacer in self.list_spacers]
        
        list_sum_list_repeats = [sum(list_length_repeats[:x]) for x in range(len(list_length_repeats) + 1)]
        list_sum_list_spacers = [sum(list_length_spacers[:x]) for x in range(len(list_length_spacers) + 1)]
        
        list_repeat_starts = [r + s for r, s in zip(list_sum_list_repeats, list_sum_list_spacers)]        
        list_repeat_ends = [r + s for r, s in zip(list_sum_list_repeats[1:], list_sum_list_spacers)]
        
        list_repeat_starts = [r_start + 1 for r_start in list_repeat_starts]
        list_repeat_ends = [r_end + 1 for r_end in list_repeat_ends]
        
        self.list_repeat_intervals = [(r_start, r_end) for 
                                      r_start, r_end in zip(list_repeat_starts, list_repeat_ends)]       
        
    def _create_crispr_seq(self):
        for r, s in zip(self.list_repeats, self.list_spacers):
            self.crispr_seq += (r + s)
        self.crispr_seq += self.list_repeats[-1]            
            
    def _create_file(self):
        with open('fasta_to_do_orf_{}.fa'.format(self.index), 'w') as f:
            f.write('>seq_to_score\n')
            f.write(self.crispr_seq)
            
    def _call_prodigal(self):
        cmd = 'tools/prodigal/prodigal'
        cmd += ' -i fasta_to_do_orf_{}.fa -o prodi_{}.txt -c -p meta'.format(self.index, self.index)

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()
        
    def _read_from_orf_file(self):                
        with open('prodi_{}.txt'.format(self.index), 'r') as f:
            lines = f.readlines()

        for index, line in enumerate(lines):
            if 'CDS' in line:
                next_line = lines[index + 1]
                region_string = re.search(r'\d+\.\.\d+', line).group(0)
                region_tuple = tuple(int(x) for x in region_string.split('..'))
                flag_complement = 1 if 'complement' in line else 0
                tuple_region_flag_complement = region_tuple[0], region_tuple[1], flag_complement
                
                match = re.search(r'conf=\d+\.\d+;', next_line)                
                confidence = float(match.group(0)[5:-1])
                
                self.dict_cds[tuple_region_flag_complement] = confidence        
        
    def _find_best_score(self):        
        crispr_length = len(self.crispr_seq)
        list_tuple_region_flag_complement = self.dict_cds.keys()        
        for tuple_region_flag_complement in list_tuple_region_flag_complement:
            if tuple_region_flag_complement[2] == 0:
                reg_start = tuple_region_flag_complement[0]
            else:
                reg_start = crispr_length - tuple_region_flag_complement[1]
            
            for repeat_interval in self.list_repeat_intervals:
                interval_start = repeat_interval[0]
                interval_end = repeat_interval[1]                
                if interval_end >= reg_start >= interval_start:
                    score = self.dict_cds[tuple_region_flag_complement]
                    if score > self.orf_best_score:
                        self.orf_best_score = score
                        
    def _clean_up(self):
        try:
            os.remove('prodi_{}.txt'.format(self.index))
            os.remove("fasta_to_do_orf_{}.fa".format(self.index))
        except OSError:
            pass        
                        
    def output(self):
        return float(self.orf_best_score)


class CrisprHmmer(RepeatSpacersFeatures):
    def __init__(self, index, list_repeats, list_spacers):
        start_time = time()

        super().__init__(index, list_repeats, list_spacers)
        self.crispr_seq = ''
        
        self.hmm_score = None
        
        self._create_crispr_seq()
        self._create_file()
        self._call_prodigal()
        self._run_hmm_search('tools/hmm_search/models_tandem.hmm')
        self._get_best_score_from_hmm()
        self._clean_up()

        end_time = time()
        time_taken = end_time - start_time
        #print(f"Hmmr score took {time_taken}")
    
    def _create_crispr_seq(self):
        for r, s in zip(self.list_repeats, self.list_spacers):
            self.crispr_seq += (r + s)
        self.crispr_seq += self.list_repeats[-1]            
            
    def _create_file(self):
        with open('fasta_to_do_hmm_{}.fa'.format(self.index), 'w') as f:
            f.write('>seq_to_score\n')
            f.write(self.crispr_seq)
            
    def _call_prodigal(self):
        cmd = 'tools/prodigal/prodigal'
        cmd += ' -i fasta_to_do_hmm_{}.fa -p meta -a protein_{}.fa'.format(self.index, self.index)

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        a, b = process.communicate()
    
    def _run_hmm_search(self, hmm_model):
        if os.stat('protein_{}.fa'.format(self.index)).st_size != 0:
            cmd = 'tools/hmm_search/hmmsearch --tblout {} {} {}'.format('result_hmm_{}.out'.format(self.index),
                                                                        hmm_model,
                                                                        'protein_{}.fa'.format(self.index))
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
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


class ConsensusBlast(ConsensusFeatures):
    def __init__(self, index, consensus_seq):
        start_time = time()

        super().__init__(index, consensus_seq)
        self.db1_score = None
        self.db2_score = None
        
        self._get_blast_scores()

        end_time = time()
        time_taken = end_time - start_time
        #print(f"Blast score took {time_taken}")
        
    def _write_file_with_consensus(self):
        with open('file_with_consensus_{}.fa'.format(self.index), 'w') as f:
            f.write('>consensus\n')
            f.write(self.consensus_seq)
            
    def _blast_vs_databeses(self):
            for number in [1, 2]:
                db_file = 'Verified_repeats_dataset{}.fa'.format(number)

                cmd = 'tools/blasting/blastn -query file_with_consensus_{}.fa'.format(self.index)
                cmd += ' -db tools/blasting/'
                cmd += db_file
                cmd += ' -word_size 6'
                cmd += ' -outfmt 6 -out output_vs_db_{}_{}.txt'.format(number, self.index)

                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                process.communicate()
    
    def _get_score_from_file(self, file_name):
            with open(file_name, 'r') as f:
                line = f.readline()
                line_numbers = line.split('\t')
                if len(line_numbers) > 0:
                    score = line.split('\t')[-1].strip()
                else:
                    score = '0.0'
                
                if score == '':
                    score = '0.0'
            return score        
        
    def _get_blast_scores(self):
        self._write_file_with_consensus()
        self._blast_vs_databeses()
        
        self.db1_score = self._get_score_from_file('output_vs_db_1_{}.txt'.format(self.index))
        self.db2_score = self._get_score_from_file('output_vs_db_2_{}.txt'.format(self.index))
        
        self._clean_up()
        
    def _clean_up(self):
        try:
            os.remove('file_with_consensus_{}.fa'.format(self.index))
            os.remove('output_vs_db_1_{}.txt'.format(self.index))
            os.remove('output_vs_db_2_{}.txt'.format(self.index))
        except OSError:
            pass          
        
    def output(self):
        return float(self.db1_score), float(self.db2_score)


class ConsensusMFEScore(ConsensusFeatures):
    def __init__(self, index, consensus_seq):

        start_time = time()

        super().__init__(index, consensus_seq)
        self.consensus_seq = consensus_seq
        
        self.mfe_score = None
        self._calculate_mfe_score()

        end_time = time()
        time_taken = end_time - start_time
        #print(f"MFE score took {time_taken}")
    
    def _calculate_mfe_score(self):        
        cmd_orig = 'echo ' + self.consensus_seq + ' | tools/rna_fold/RNAfold --noLP > out_st_{}.txt'.format(self.index)
        process = subprocess.Popen(cmd_orig, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        with open('out_st_{}.txt'.format(self.index), 'r') as f:
            list_lines = f.readlines()
            try:
                last_line = list_lines[-1]
                search = re.search(r'-*\d+.\d+', last_line)
                self.mfe_score = search.group(0)
            except IndexError:
                self.mfe_score = 0.0

        try:
            os.remove('out_st_{}.txt'.format(self.index))
        except OSError:
            pass

        try:
            os.remove("rna.ps")
        except OSError:
            pass
    
    def output(self):
        return float(self.mfe_score)


class ConsensusMFEScoreNew(ConsensusFeatures):
    def __init__(self, index, consensus_seq):

        start_time = time()

        super().__init__(index, consensus_seq)
        self.consensus_seq = consensus_seq

        self.mfe_score = None
        self._calculate_mfe_score()

        end_time = time()
        time_taken = end_time - start_time
        #print(f"MFE score took {time_taken}")

    def _calculate_mfe_score(self):
        cmd_orig = 'echo ' + self.consensus_seq + ' | tools/rna_fold/RNAfold --noLP --noPS'.format(self.index)
        process = subprocess.Popen(cmd_orig, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        result, errors = process.communicate()

        result = str(result)
        result = result.splitlines()

        try:
            last_line = result[-1]
            search = re.search(r'-*\d+.\d+', last_line)
            self.mfe_score = search.group(0)
        except IndexError:
            self.mfe_score = 0.0

    def output(self):
        return float(self.mfe_score)


class EdenScore:
    def __init__(self, consensus_sequence, eden_classifier):
        self.consensus_sequence = consensus_sequence
        self.classifier = eden_classifier
        
    def output(self):
        return self.classifier.score_input([("seq", self.consensus_sequence)])[0]


class FeatureExtractor:
    def __init__(self, index, crispr_candidate, eden_classifier, features_to_extract=None):
        self.index = index
        self.crispr_candidate = crispr_candidate        
        self.eden_classifier = eden_classifier
        self.features_to_extract = features_to_extract
        self.list_repeats = crispr_candidate.list_repeats
        self.list_repeats_gaped = crispr_candidate.list_repeats_gaped
        self.list_spacers = crispr_candidate.list_spacers
        
        self.dict_features = {}
        
    def extract(self):
        consensus = self.crispr_candidate.consensus        
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
            repeat_similarity, spacer_similarity = CrisprSimilarityNew(self.index, self.list_repeats, self.list_spacers).output()
            self.dict_features['repeat_similarity'] = repeat_similarity
            self.dict_features['spacer_similarity'] = spacer_similarity
            
        if 'orf_score' in self.features_to_extract:
            orf_score = CrisprORFScore(self.index, self.list_repeats, self.list_spacers).output()
            self.dict_features['orf_score'] = orf_score
            
        if 'hmmr_score' in self.features_to_extract:
            hmmr_score = CrisprHmmer(self.index, self.list_repeats, self.list_spacers).output()
            self.dict_features['hmmr_score'] = hmmr_score
            
        if ('blast_score_1' in self.features_to_extract) or ('blast_score_2' in self.features_to_extract):
            blast_score_1, blast_score_2 = ConsensusBlast(self.index, consensus).output()
            self.dict_features['blast_score_1'] = blast_score_1
            self.dict_features['blast_score_2'] = blast_score_2
            
        if 'mfe_score' in self.features_to_extract:
            mfe_score = ConsensusMFEScoreNew(self.index, consensus).output()
            self.dict_features['mfe_score'] = mfe_score
            
        if 'eden_score' in self.features_to_extract:
            eden_score = EdenScore(consensus, self.eden_classifier).output()
            self.dict_features['eden_score'] = eden_score
        
        feature_vector = [self.dict_features[feature] for feature in self.features_to_extract]
        
        return np.asarray(feature_vector).reshape(1, -1)
 
