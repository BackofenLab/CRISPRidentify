import numpy as np

from components.components_evaluation import BulkFeatureExtractor
from components.components_evaluation import FeatureExtractor
from components.components_evaluation import get_full_vector
from components.components_detection_refinement import AdvancedFuzzySearchFilter


class ArrayEvaluation:
    def __init__(self, dict_crispr_array_candidates, list_ml_classifiers, list_features, parameters, flag_dev_mode):
        self.dict_crispr_array_candidates = dict_crispr_array_candidates
        self.list_ml_classifiers = list_ml_classifiers
        self.list_features = list_features
        self.parameters = parameters
        self.flag_dev_mode = flag_dev_mode

        self.dict_scored_result = {}
        self.dict_scored_result_with_all_vectors = {}

        self.dict_bona_fide = {}
        self.dict_alternative = {}
        self.dict_possible = {}
        self.dict_possible_discarded = {}
        self.dict_low_score = {}

        self._load_filter()
        self._extract_features_and_evaluate()
        self._split_into_categories()

    def _load_filter(self):
        self.param_min_avg_repeat_length = self.parameters["param_min_avg_repeat_length"]
        self.param_max_avg_repeat_length = self.parameters["param_max_avg_repeat_length"]
        self.param_min_avg_spacer_length = self.parameters["param_min_avg_spacer_length"]
        self.param_max_avg_spacer_length = self.parameters["param_max_avg_spacer_length"]
        self.param_min_repeats = self.parameters["param_min_repeats"]
        self.param_max_identical_spacers = self.parameters["param_max_identical_spacers"]
        self.param_max_identical_cluster_spacers = self.parameters["param_max_identical_cluster_spacers"]
        self. afsf = AdvancedFuzzySearchFilter(min_column_dominance_repeat=0.6,
                                               max_spacer_length=140, max_column_dominance_spacer=0.8,
                                               max_allowed_consecutive_spacers=self.param_max_identical_cluster_spacers,
                                               max_allowed_same_spacers=self.param_max_identical_spacers,
                                               max_inconsistent_columns=5,
                                               min_avg_repeat_length=self.param_min_avg_repeat_length,
                                               max_avg_repeat_length=self.param_max_avg_repeat_length,
                                               min_avg_spacer_length=self.param_min_avg_spacer_length,
                                               max_avg_spacer_length=self.param_max_avg_spacer_length,
                                               min_repeats=self.param_min_repeats)

    def _extract_features_and_evaluate(self):
        bfe = BulkFeatureExtractor(self.dict_crispr_array_candidates)
        results = bfe.output()
        blast_results, orf_results, hmm_results, mfe_results = results
        blast_scores_1, blast_scores_2 = blast_results

        list_features = ['repeat_len', 'number_repeats', 'repeat_similarity',
                         'at_richness', 'avg_spacer_len', 'spacer_similarity',
                         'number_mismatches', 'spacer_evenness']

        for key, list_crispr_candidates in self.dict_crispr_array_candidates.items():
            self.dict_scored_result[key] = []
            self.dict_scored_result_with_all_vectors[key] = []
            for index, crispr_candidate in enumerate(list_crispr_candidates):
                final_score = 0

                feature_vector = FeatureExtractor(0, crispr_candidate, list_features).extract()[0]

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

                final_score = final_score / len(self.list_ml_classifiers)
                score_crispr_candidate_feature_list = [final_score, crispr_candidate, feature_vectors]
                self.dict_scored_result[key].append(score_crispr_candidate_feature_list)

                all_feature_vectors = [feature_vector_8, feature_vector_9, feature_vector_10]
                score_crispr_candidate_all_feature_tuple = final_score, crispr_candidate, all_feature_vectors
                self.dict_scored_result_with_all_vectors[key].append(score_crispr_candidate_all_feature_tuple)

    def _split_into_categories(self):
        for key, data in self.dict_scored_result.items():
            data_pre_possible = [candidate for candidate in data if 0.75 > candidate[0] >= 0.5]
            data_alternative = [candidate for candidate in data if candidate[0] >= 0.75]
            data_alternative_filtered = []
            data_bad = [candidate for candidate in data if candidate[0] < 0.5]

            if data_alternative:
                for element in data_alternative:
                    crispr = element[1]
                    if self.afsf(crispr):
                        data_alternative_filtered.append(element)
                    else:
                        data_bad.append(element)

                if data_alternative_filtered:
                    data_alternative_filtered = sorted(data_alternative_filtered, key=lambda x: x[0], reverse=True)
                    best_candidate = data_alternative_filtered[0]
                    data_alternative_filtered = data_alternative_filtered[1:]

                    self.dict_bona_fide[key] = [best_candidate]
                    if data_alternative_filtered:
                        self.dict_alternative[key] = data_alternative_filtered

            if data_pre_possible:
                if key in self.dict_bona_fide:
                    data_show_in_alternative = [candidate for candidate in data_pre_possible if candidate[0] >= 0.6]
                    if data_show_in_alternative:
                        data_show_in_alternative_filtered = []
                        for element in data_show_in_alternative:
                            crispr = element[1]
                            if self.afsf(crispr):
                                data_show_in_alternative_filtered.append(element)
                            else:
                                data_bad.append(element)

                        if key in self.dict_alternative:
                            self.dict_alternative[key] += data_show_in_alternative_filtered
                        else:
                            self.dict_alternative[key] = data_show_in_alternative_filtered

                else:
                    data_pre_possible = sorted(data_pre_possible, key=lambda x: x[0], reverse=True)
                    best_possible_candidate = data_pre_possible[0]
                    possible_discarded = data_pre_possible[1:]

                    if self.afsf(best_possible_candidate[1]):
                        self.dict_possible[key] = [best_possible_candidate]
                    else:
                        data_bad.append(best_possible_candidate)

                    if possible_discarded:
                        self.dict_possible_discarded[key] = possible_discarded

            if data_bad:
                self.dict_low_score[key] = data_bad

    def _split_into_categories_with_additional_classifier(self):

        for key, data in self.dict_scored_result_with_all_vectors.items():
            data_pre_possible = [candidate for candidate in data if 0.75 > candidate[0] >= 0.5]
            data_alternative = [candidate for candidate in data if candidate[0] >= 0.75]
            data_alternative_filtered = []
            data_bad = [candidate for candidate in data if candidate[0] < 0.5]

            if data_bad:
                self.dict_low_score[key] = data_bad

            if self.flag_possible_differential_model == "possible":
                if data_alternative:
                    data_alternative = sorted(data_alternative, key=lambda x: x[0], reverse=True)
                    best_candidate = data_alternative[0]
                    data_alternative = data_alternative[1:]

                    self.dict_bona_fide[key] = best_candidate
                    if data_alternative:
                        self.dict_alternative[key] = data_alternative
            else:
                if data_alternative:
                    for element in data_alternative:
                        crispr = element[1]
                        if self.afsf(crispr):
                            data_alternative_filtered.append(element)
                        else:
                            data_pre_possible.append(element)

                        data_alternative_filtered = sorted(data_alternative_filtered, key=lambda x: x[0], reverse=True)
                        best_candidate = data_alternative_filtered[0]
                        data_alternative_filtered = data_alternative_filtered[1:]

                        self.dict_bona_fide[key] = [best_candidate]
                        if data_alternative_filtered:
                            self.dict_alternative[key] = data_alternative_filtered
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
                        self.dict_bona_fide[key] = [best_candidate]
                        alternative = data_alternative_sorted[1:]
                        if alternative:
                            self.dict_alternative[key] = alternative
                    else:
                        self.dict_bona_fide[key] = [best_candidate_prev_model]
                        if data_alternative_prev_model:
                            self.dict_alternative[key] = data_alternative_prev_model

            if data_pre_possible:
                data_pre_possible = sorted(data_pre_possible, key=lambda x: x[0], reverse=True)

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
                    possible_discarded = data_pre_possible_sorted
                    self.dict_possible_discarded[key] = possible_discarded

    def output(self):
        return [self.dict_bona_fide, self.dict_alternative, self.dict_possible,
                self.dict_possible_discarded, self.dict_low_score]

