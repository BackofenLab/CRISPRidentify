import pickle
import os
import sys
from os.path import basename

from detection import Detector
from refine_steps import Refine
from feature_extraction import FeatureExtractor
from bulk_feature_extraction import BulkFeatureExtractor
from feature_extraction import CrisprNumberMismatches
from non_crispr_computations import AdditionalCalculations
from array_refinder_insertions_deletions import ArrayRefinerInsertionsDeletions
from array_refiner_isolated_repeat import ArrayRefinerIsolatedRepeat
from summary_maker import SummaryMaker
import multiprocessing
from multiprocessing import Pool

from components_result_parser import ResultParser
from to_complete_vector import get_full_vector
from time import time
import numpy as np


class Predictor(object):
    def __init__(self, result_folder_path, file_path, eden_classifier, list_ml_classifiers,
                 possible_differentiate_model, flag_possible_differential_model,
                 list_features, flag_is, flag_cas, flag_degenerated, flag_parallel,
                 flag_cpu, flag_fast_run, flag_enhancement_max_min, flag_enhancement_start_end,
                 parameters, log_file):
        self.parameters = parameters
        self.flag_parallel = flag_parallel
        self.flag_cpu = flag_cpu
        self.result_folder_path = result_folder_path
        self.file_path = file_path
        self.eden_classifier = eden_classifier
        self.list_ml_classifiers = list_ml_classifiers
        self.possible_differentiate_model = possible_differentiate_model
        self.flag_possible_differential_model = flag_possible_differential_model
        self.list_features = [features.strip().split(".") for features in list_features]
        self.flag_is = flag_is
        self.flag_cas = flag_cas
        self.flag_degenerated = flag_degenerated
        self.flag_fast_run = flag_fast_run
        self.flag_enhancement_max_min = flag_enhancement_max_min
        self.flag_enhancement_start_end = flag_enhancement_start_end
        self.log_file = log_file
        
        self.dict_detection_output = {}
        self.dict_filtered_detection = {}
        self.dict_refined_detection = {}
        self.dict_scored_result = {}
        self.dict_scored_result_with_all_vectors = {}
        self.dict_crispr_indexes = {}
        
        self._detect_data()
        self._filter_detection()
        self._refine_detection()
        self._refine_isolated_repeats()
        if self.flag_parallel:
            #self._score_data_parallel()
            self._score_data()  # Temporary fix since the bulk feature extraction is used
        else:
            self._score_data()
        if flag_possible_differential_model is None:
            self._create_categories()
        else:
            self._create_categories_with_additional_classifier()

        self._refine_insertions_deletions_regex()
        self._index_crispr_candidates()
        self._write_output()
        self._write_summary()
        self._write_csv_summary()
        self._write_bed_summary()
        
    def _detect_data(self):
        start_time = time()
        self.detector = Detector(self.file_path, self.flag_parallel, self.flag_cpu,
                                 self.flag_fast_run, self.flag_enhancement_max_min, self.flag_enhancement_start_end,
                                 self.parameters, self.log_file)

        self.dict_detection_output = self.detector.output()
        end_time = time()
        time_taken = end_time - start_time
        #print(f"Detection step took {time_taken}")
        
    def _filter_detection(self):
        """
        Filtering cases there is a lot of mismatches
        And cases where where are no repeats or spacers
        """
        start_time = time()

        for key, list_candidates in self.dict_detection_output.items():
            list_filtered_candidates = []
            for index, candidate in enumerate(list_candidates):
                list_gaped_repeats = candidate.list_repeats_gaped
                list_spacers = candidate.list_spacers
                if any(list_gaped_repeats) and any(list_spacers):
                    consensus_gaped = candidate.consensus_gaped
                    mismatches = CrisprNumberMismatches(index, list_gaped_repeats, consensus_gaped,).output()
                    if mismatches/len(list_gaped_repeats) < 3:
                        list_filtered_candidates.append(candidate)
            
            self.dict_filtered_detection[key] = list_filtered_candidates

        end_time = time()
        time_taken = end_time - start_time
        #print(f"Filtering detection step took {time_taken}")

    def _refine_detection(self):

        start_time = time()

        refine = Refine()
        for key, list_candidates in self.dict_filtered_detection.items():
            list_additional_candidates = refine.search_for_part_spacers_repeat(list_candidates, 6)
            list_candidates += list_additional_candidates
            self.dict_refined_detection[key] = list_candidates

        end_time = time()
        time_taken = end_time - start_time
        #print(f"Refine detection step took {time_taken}")

    def _refine_isolated_repeats(self):
        dict_refined_isolated_repeats = {}

        for key, list_candidates in self.dict_refined_detection.items():
            list_refined_isolated = [ArrayRefinerIsolatedRepeat(cand, 140).output() for cand in list_candidates]
            dict_refined_isolated_repeats[key] = list_refined_isolated

        self.dict_refined_detection = dict_refined_isolated_repeats

    def _score_data(self):
        print('Extracting features for {} candidates'.format(sum(len(self.dict_refined_detection[key])
                                                                 for key in self.dict_refined_detection)))

        bfe = BulkFeatureExtractor(self.dict_refined_detection)
        results = bfe.output()
        blast_results, orf_results, hmm_results, mfe_results = results
        blast_scores_1, blast_scores_2 = blast_results

        list_features = ['repeat_len', 'number_repeats', 'repeat_similarity',
                         'at_richness', 'avg_spacer_len', 'spacer_similarity',
                         'number_mismatches', 'spacer_evenness']

        for key, list_crispr_candidates in self.dict_refined_detection.items():
            self.dict_scored_result[key] = []
            self.dict_scored_result_with_all_vectors[key] = []
            for index, crispr_candidate in enumerate(list_crispr_candidates):
                final_score = 0

                feature_vector = FeatureExtractor(0, crispr_candidate, self.eden_classifier,
                                                  list_features).extract()[0]

                mfe = mfe_results[key][index]
                orf = orf_results[key][index]
                hmmr = hmm_results[key][index]
                blast1 = blast_scores_1[key][index]
                blast2 = blast_scores_2[key][index]

                feature_vector_8_incomplete = feature_vector[np.array([2, 4, 5, 6, 7])]
                rest_8 = np.asarray([mfe, orf, blast1])
                feature_vector_8 = np.concatenate((feature_vector_8_incomplete, rest_8))
                feature_vector_8 = feature_vector_8.reshape(1, -1)

                feature_vector_9_incomplete = feature_vector[np.array([1, 2, 4, 5, 7])]
                rest_9 = np.asarray([mfe, orf, hmmr, blast2])
                feature_vector_9 = np.concatenate((feature_vector_9_incomplete, rest_9))
                feature_vector_9 = feature_vector_9.reshape(1, -1)

                feature_vector_10_incomplete = feature_vector[np.array([0, 2, 3, 4, 5, 6, 7])]
                rest_10 = np.asarray([hmmr, blast1, blast2])
                feature_vector_10 = np.concatenate((feature_vector_10_incomplete, rest_10))
                feature_vector_10 = feature_vector_10.reshape(1, -1)

                dict_feature_vectors = {8: feature_vector_8,
                                        9: feature_vector_9,
                                        10: feature_vector_10}

                feature_vectors = []
                for ml_classifier, feature_names in zip(self.list_ml_classifiers, self.list_features):
                    len_features = len(feature_names)
                    feature_vector = dict_feature_vectors[len_features]
                    feature_vectors.append(feature_vector)
                    final_score += ml_classifier.predict_proba(feature_vector)[0][1]

                final_score = final_score/len(self.list_ml_classifiers)
                score_crispr_candidate_feature_tuple = final_score, crispr_candidate, feature_vectors
                self.dict_scored_result[key].append(score_crispr_candidate_feature_tuple)

                all_feature_vectors = [feature_vector_8, feature_vector_9, feature_vector_10]
                score_crispr_candidate_all_feature_tuple = final_score, crispr_candidate, all_feature_vectors
                self.dict_scored_result_with_all_vectors[key].append(score_crispr_candidate_all_feature_tuple)

    @staticmethod
    def _parallel_feature_extraction(input_for_feature_extraction):
        index, crispr_candidate, feature_names = input_for_feature_extraction
        return FeatureExtractor(index, crispr_candidate, None, feature_names).extract()
    
    def _create_categories(self):
        self.dict_best = {}
        self.dict_alternative = {}
        self.dict_possible = {}
        self.dict_possible_discarded = {}
        self.dict_bad = {}
        
        for key, data in self.dict_scored_result.items():
            data_pre_possible = [candidate for candidate in data if 0.75 > candidate[0] >= 0.5]
            data_alternative = [candidate for candidate in data if candidate[0] >= 0.75]
            data_bad = [candidate for candidate in data if candidate[0] < 0.5]
                
            if data_bad:
                self.dict_bad[key] = data_bad
                
            if data_alternative:
                data_alternative = sorted(data_alternative, key=lambda x: x[0], reverse=True)
                best_candidate = data_alternative[0]
                data_alternative = data_alternative[1:]
                
                self.dict_best[key] = best_candidate
                if data_alternative:
                    self.dict_alternative[key] = data_alternative

            if data_pre_possible:
                data_pre_possible = sorted(data_pre_possible, key=lambda x: x[0], reverse=True)
                best_possible_candidate = data_pre_possible[0]
                possible_discarded = data_pre_possible[1:]
                self.dict_possible[key] = [best_possible_candidate]
                self.dict_possible_discarded[key] = possible_discarded

    def _create_categories_with_additional_classifier(self):
        self.dict_best = {}
        self.dict_alternative = {}
        self.dict_possible = {}
        self.dict_possible_discarded = {}
        self.dict_bad = {}

        for key, data in self.dict_scored_result_with_all_vectors.items():
            data_pre_possible = [candidate for candidate in data if 0.75 > candidate[0] >= 0.5]
            data_alternative = [candidate for candidate in data if candidate[0] >= 0.75]
            data_bad = [candidate for candidate in data if candidate[0] < 0.5]

            if data_bad:
                self.dict_bad[key] = data_bad

            if self.flag_possible_differential_model == "possible":
                if data_alternative:
                    data_alternative = sorted(data_alternative, key=lambda x: x[0], reverse=True)
                    best_candidate = data_alternative[0]
                    data_alternative = data_alternative[1:]

                    self.dict_best[key] = best_candidate
                    if data_alternative:
                        self.dict_alternative[key] = data_alternative
            else:
                if data_alternative:
                    data_alternative = sorted(data_alternative, key=lambda x: x[0], reverse=True)
                    best_candidate_prev_model = data_alternative[0]
                    data_alternative_prev_model = data_alternative[1:]

                    vectors_alternative = [get_full_vector(data[2]) for data in data_alternative]
                    scores_new_model = [self.possible_differentiate_model.predict_proba(v)[0][1] for v in
                                        vectors_alternative]

                    scores_new_model, data_alternative_sorted = zip(*sorted(zip(scores_new_model, data_alternative),
                                                                    key=lambda x: x[0], reverse=True))

                    best_candidate = data_alternative_sorted[0]
                    best_score = scores_new_model[0]
                    label = 1.0 if best_score >= 0.5 else 0.0

                    if label == 1.0:
                        self.dict_best[key] = best_candidate
                        alternative = data_alternative_sorted[1:]
                        if alternative:
                            self.dict_alternative[key] = alternative
                    else:
                        self.dict_best[key] = best_candidate_prev_model
                        if data_alternative_prev_model:
                            self.dict_alternative[key] = data_alternative_prev_model

            if data_pre_possible:
                data_pre_possible = sorted(data_pre_possible, key=lambda x: x[0], reverse=True)

                #best_possible_candidate_prev_model = data_pre_possible[0]
                #possible_discarded_prev_model = data_pre_possible[1:]

                vectors_pre_possible = [get_full_vector(data[2]) for data in data_pre_possible]
                scores_new_model = [self.possible_differentiate_model.predict_proba(v)[0][1] for v in vectors_pre_possible]

                scores_new_model, data_pre_possible_sorted = zip(*sorted(zip(scores_new_model, data_pre_possible),
                                                                 key=lambda x: x[0], reverse=True))

                best_possible_candidate = data_pre_possible_sorted[0]
                best_score = scores_new_model[0]
                label = 1.0 if best_score >= 0.5 else 0.0

                if label == 1.0:
                    self.dict_possible[key] = [best_possible_candidate]
                    possible_discarded = data_pre_possible_sorted[1:]
                    self.dict_possible_discarded[key] = possible_discarded
                else:
                    #self.dict_possible[key] = [best_possible_candidate_prev_model]
                    possible_discarded = data_pre_possible_sorted
                    self.dict_possible_discarded[key] = possible_discarded

    def _refine_insertions_deletions_regex(self):
        for key, value in self.dict_best.items():
            val_copy = list(value[:])
            crispr_candidate = value[1]
            ref_candidate = ArrayRefinerInsertionsDeletions(crispr_candidate).output()
            val_copy[1] = ref_candidate
            val_copy = tuple(val_copy)
            self.dict_best[key] = val_copy

    def _index_crispr_candidates(self):
        """
        for index, key in enumerate(sorted(self.dict_best.keys()), 1):
            self.dict_crispr_indexes[key] = index

        for index, key in enumerate(sorted(self.dict_possible.keys()), len(self.dict_best) + 1):
            if key not in self.dict_crispr_indexes:
                self.dict_crispr_indexes[key] = index

        for index, key in enumerate(sorted(self.dict_bad.keys()),
                                    len(self.dict_best) + len(self.dict_possible) + 1):
            if key not in self.dict_crispr_indexes:
                self.dict_crispr_indexes[key] = index

        """

        for index, key in enumerate(sorted(self.dict_best.keys()), 1):
            self.dict_crispr_indexes[key] = index

        index = len(self.dict_best)
        cur_index = index + 1

        for key in sorted(self.dict_possible.keys()):
            if key not in self.dict_crispr_indexes:
                self.dict_crispr_indexes[key] = cur_index
                cur_index += 1

        for key in sorted(self.dict_possible_discarded.keys()):
            if key not in self.dict_crispr_indexes:
                self.dict_crispr_indexes[key] = cur_index
                cur_index += 1

        for key in sorted(self.dict_bad.keys()):
            if key not in self.dict_crispr_indexes:
                self.dict_crispr_indexes[key] = cur_index
                cur_index += 1

    def _write_output(self):
        self.file_base = basename(self.file_path)
                
        if not os.path.exists(self.result_folder_path):
            os.makedirs(self.result_folder_path)
            
        result_path = self.result_folder_path + "/Result_" + self.file_base.split(".")[0]
        if not os.path.exists(result_path):
            os.makedirs(result_path)       
        
        print("Writing Crispr Candidates into files")
        f_best = open(result_path + '/Bona_fide_Candidates.txt', 'w')
        f_alternative = open(result_path + '/Alternative_Candidates.txt', 'w')
        f_possible = open(result_path + '/Possible_Candidates.txt', 'w')
        f_possible_discarded = open(result_path + '/Possible_Discarded_Candidates.txt', 'w')
        f_bad = open(result_path + '/Low_score_Candidates.txt', 'w')
        
        for key in sorted(self.dict_best.keys()):
            crispr = self.dict_best[key][1]
            crispr_stats = crispr.compute_stats()
            f_best.write("CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                         .format(self.dict_crispr_indexes[key], crispr_stats["start"], crispr_stats["end"],
                                 crispr_stats["number_repeats"], crispr_stats["avg_repeat"],
                                 crispr_stats["avg_spacer"]))

            f_best.write(crispr.dot_repr())

            f_best.write("\n")
            list_reported_features = []
            for index, feature_list in enumerate(self.list_features):
                for feature, value in zip(feature_list, self.dict_best[key][2][index][0]):
                    if feature not in list_reported_features:
                        f_best.write("{}: {}\n".format(feature, value))
                        list_reported_features.append(feature)

            f_best.write("\n")
            f_best.write("Certainty Score: {}\n\n".format(self.dict_best[key][0]))
            f_best.write('\n{}\n\n'.format('=' * 100))
            
        for key in sorted(self.dict_alternative.keys()):
            for candidate in self.dict_alternative[key]:
                score = candidate[0]
                crispr = candidate[1]
                crispr_stats = crispr.compute_stats()
                f_alternative.write("Alternative CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                                    .format(self.dict_crispr_indexes[key], crispr_stats["start"], crispr_stats["end"],
                                            crispr_stats["number_repeats"], crispr_stats["avg_repeat"],
                                            crispr_stats["avg_spacer"]))

                f_alternative.write(crispr.dot_repr())

                f_alternative.write("\n")
                list_reported_features = []
                for index, feature_list in enumerate(self.list_features):
                    for feature, value in zip(feature_list, candidate[2][index][0]):
                        if feature not in list_reported_features:
                            f_alternative.write("{}: {}\n".format(feature, value))
                            list_reported_features.append(feature)

                f_alternative.write("\n")
                f_alternative.write("Certainty Score: {}\n\n\n\n".format(score))

            f_alternative.write('\n{}\n\n'.format('=' * 100))

        for key in sorted(self.dict_possible.keys()):
            for candidate in self.dict_possible[key]:
                score = candidate[0]
                crispr = candidate[1]
                feature_vector = candidate[2][0]
                crispr_stats = crispr.compute_stats()
                f_possible.write("Possible CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                                 .format(self.dict_crispr_indexes[key], crispr_stats["start"], crispr_stats["end"],
                                         crispr_stats["number_repeats"], crispr_stats["avg_repeat"],
                                         crispr_stats["avg_spacer"]))
                f_possible.write(crispr.dot_repr())

                f_possible.write("\n")
                list_reported_features = []
                for index, feature_list in enumerate(self.list_features):
                    for feature, value in zip(feature_list, candidate[2][index][0]):
                        if feature not in list_reported_features:
                            f_possible.write("{}: {}\n".format(feature, value))
                            list_reported_features.append(feature)

                f_possible.write("\n")
                f_possible.write("Certainty Score: {}\n\n\n\n".format(score))

            f_possible.write('\n{}\n\n'.format('=' * 100))

        for key in sorted(self.dict_possible_discarded.keys()):
            for candidate in self.dict_possible_discarded[key]:
                score = candidate[0]
                crispr = candidate[1]
                feature_vector = candidate[2][0]
                crispr_stats = crispr.compute_stats()
                f_possible_discarded.write("Possible Discarded CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                                 .format(self.dict_crispr_indexes[key], crispr_stats["start"], crispr_stats["end"],
                                         crispr_stats["number_repeats"], crispr_stats["avg_repeat"],
                                         crispr_stats["avg_spacer"]))
                f_possible_discarded.write(crispr.dot_repr())

                f_possible_discarded.write("\n")
                list_reported_features = []
                for index, feature_list in enumerate(self.list_features):
                    for feature, value in zip(feature_list, candidate[2][index][0]):
                        if feature not in list_reported_features:
                            f_possible_discarded.write("{}: {}\n".format(feature, value))
                            list_reported_features.append(feature)

                f_possible_discarded.write("\n")
                f_possible_discarded.write("Certainty Score: {}\n\n\n\n".format(score))

            f_possible_discarded.write('\n{}\n\n'.format('=' * 100))
        
        for key in sorted(self.dict_bad.keys()):
            for candidate in self.dict_bad[key]:
                score = candidate[0]
                crispr = candidate[1]
                crispr_stats = crispr.compute_stats()
                f_bad.write("Low_score CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                            .format(self.dict_crispr_indexes[key], crispr_stats["start"], crispr_stats["end"],
                                    crispr_stats["number_repeats"], crispr_stats["avg_repeat"],
                                    crispr_stats["avg_spacer"]))
                f_bad.write(crispr.dot_repr())

                f_bad.write("\n")
                list_reported_features = []
                for index, feature_list in enumerate(self.list_features):
                    for feature, value in zip(feature_list, candidate[2][index][0]):
                        if feature not in list_reported_features:
                            f_bad.write("{}: {}\n".format(feature, value))
                            list_reported_features.append(feature)

                f_bad.write("\n")
                f_bad.write("Certainty Score: {}\n\n\n\n".format(score))

            f_bad.write('\n{}\n\n'.format('=' * 100))
        
        f_best.close()
        f_possible.close()
        f_possible_discarded.close()
        f_alternative.close()
        f_bad.close()

    def _write_summary(self):
        acc_num = self.file_base.split(".")[0]
        result_path = self.result_folder_path + "/Result_" + self.file_base.split(".")[0]
        list_best_crisprs = [tuple_info[1] for tuple_info in self.dict_best.values()]
        list_scores = [tuple_info[0] for tuple_info in self.dict_best.values()]
        list_feature_vectors = [tuple_info[2] for tuple_info in self.dict_best.values()]
        list_indexes = [x for x in range(len(list_feature_vectors))]

        if list_best_crisprs:
            list_best_crisprs, list_indexes = (list(t) for t in zip(*sorted(zip(list_best_crisprs, list_indexes),
                                                                    key=lambda x: x[0].list_repeat_starts[0])))

            list_indexes, list_scores = (list(t) for t in zip(*sorted(zip(list_indexes, list_scores))))
            list_indexes, list_feature_vectors = (list(t) for t in zip(*sorted(zip(list_indexes,
                                                                                   list_feature_vectors))))

        full_dna = self.detector.dna
        input_file = self.file_path
        print(input_file)
        with open(input_file, "r") as f:
            header = f.readline()

        hmm_model_is_element = "tools/hmm_search/models_is_element.hmm"
        hmm_model_cas_genes = "tools/hmm_search/models_cas_genes.hmm"

        ac = AdditionalCalculations(input_file=input_file, full_dna=full_dna, list_of_crisprs=list_best_crisprs,
                                    hmm_model_is_element=hmm_model_is_element, hmm_model_cas_genes=hmm_model_cas_genes)


        summary_maker = SummaryMaker(header=header, list_crisprs=list_best_crisprs, list_scores=list_scores,
                                     list_feature_vectors=list_feature_vectors, list_feature_names=self.list_features,
                                     additional_calculations=ac, flag_is=self.flag_is, flag_cas=self.flag_cas,
                                     flag_degenerated=self.flag_degenerated, acc_num=acc_num)


        summary_maker.write_text_summary(result_path + '/Summary.txt')
        #summary_maker.write_bed_summary(result_path + '/Bed_summary.txt')
        #summary_maker.write_merged_summary(result_path + '/Merged_Summary.txt')

    def _write_csv_summary(self):
        result_path = self.result_folder_path + "/Result_" + self.file_base.split(".")[0]
        bad_path = result_path + "/Low_score_Candidates.txt"
        possible_path = result_path + "/Possible_Candidates.txt"
        alternative_path = result_path + "/Alternative_Candidates.txt"
        best_path = result_path + "/Bona_fide_Candidates.txt"

        results_best = ResultParser(best_path).output()
        results_alternative = ResultParser(alternative_path).output()
        results_possible = ResultParser(possible_path).output()
        results_bad = ResultParser(bad_path).output()

        result_csv_path = result_path + "/csv_summary.csv"
        with open(result_csv_path, "w") as f:
            f.write(",".join(["ID", "Start", "End", "Length", "Consensus repeat", "Repeat Length",
                              "Average Spacer Length", "Number of spacers", "Strand", "Category"]))
            f.write("\n")

            for result in [results_best, results_alternative, results_possible, results_bad]:
                for crispr_array_parsing_container in result:
                    category = crispr_array_parsing_container.category
                    crispr_index = crispr_array_parsing_container.crispr_index
                    start = crispr_array_parsing_container.start
                    end = crispr_array_parsing_container.end
                    number_of_repeats = crispr_array_parsing_container.number_of_repeats
                    avg_spacer_length = crispr_array_parsing_container.avg_spacer_length
                    consensus_repeat = crispr_array_parsing_container.consensus_repeat
                    strand = crispr_array_parsing_container.strand

                    string_to_write = ",".join([crispr_index, start, end, str(abs(int(end) - int(start) + 1)),
                                                consensus_repeat, str(len(consensus_repeat)), avg_spacer_length,
                                                str(int(number_of_repeats) - 1), strand, category])

                    f.write(string_to_write)
                    f.write("\n")

    def _write_bed_summary(self):
        result_path = self.result_folder_path + "/Result_" + self.file_base.split(".")[0]
        bad_path = result_path + "/Low_score_Candidates.txt"
        possible_path = result_path + "/Possible_Candidates.txt"
        alternative_path = result_path + "/Alternative_Candidates.txt"
        best_path = result_path + "/Bona_fide_Candidates.txt"

        results_best = ResultParser(best_path).output()
        results_alternative = ResultParser(alternative_path).output()
        results_possible = ResultParser(possible_path).output()
        results_bad = ResultParser(bad_path).output()

        for result, file_name in zip([results_best, results_alternative, results_possible, results_bad],
                                     ["Bona_fide_Candidates.txt", "Alternative_candidates.txt",
                                      "Possible_candidates.txt", "Low_score_Candidates.txt"]):

            bed_file_path = result_path + "/Bed_" + file_name
            with open(bed_file_path, "w") as f:
                for crispr_array_parsing_container in result:
                    acc_num = self.file_base.split(".")[0]
                    category = crispr_array_parsing_container.category
                    crispr_index = crispr_array_parsing_container.crispr_index
                    start = crispr_array_parsing_container.start
                    end = crispr_array_parsing_container.end
                    number_of_repeats = crispr_array_parsing_container.number_of_repeats
                    avg_spacer_length = crispr_array_parsing_container.avg_spacer_length
                    consensus_repeat = crispr_array_parsing_container.consensus_repeat
                    strand = crispr_array_parsing_container.strand
                    string_to_write = "\t".join([acc_num, start, end, "CRISPR" + crispr_index, ".", strand])
                    f.write(string_to_write)
                    f.write("\n")
                    
    def report_best(self):
        dict_best = {}
        
        for key, data in self.dict_scored_result.items():            
            data = sorted(data, key=lambda x: x[0])
            best_candidate = data[0:1]
            if best_candidate:                
                best_candidate_score = best_candidate[0][0]
                if best_candidate_score >= 0.75:
                    dict_best[key] = best_candidate[0][1]
        
        return dict_best
   
    def report_into_file(self, file_name):
        """Function to keep track of all data produced with one run over
        the input dataset"""
        try:
            dict_data = pickle.load(open(file_name, "rb"))
        except IOError:
            dict_data = {}
        
        key = self.file_base
        dict_data[key] = [self.dict_best, self.dict_alternative, self.dict_possible,
                          self.dict_possible_discarded, self.dict_bad]
        
        pickle.dump(dict_data, open(file_name, "wb"))
        
    def report_into_separate_file(self, folder_name):
        try:
            os.mkdir(folder_name)
        except OSError:
                pass
        key = self.file_base
        dict_data = {'best': self.dict_best,
                     'alternative': self.dict_alternative,
                     'possible': self.dict_possible,
                     'possible_discarded': self.dict_possible_discarded,
                     'bad': self.dict_bad}
        pickle.dump(dict_data, open(folder_name + '/' + key + '.pkl', "wb"))
        print("Done writing the candidates")
            
        
        
        
                
            
        

