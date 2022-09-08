from components.components_detection_refinement import SameStartEndFilter
from components.components_detection_refinement import AdvancedFuzzySearchFilter
from components.components_detection_refinement import CrisprCandidate


class DetectionRefinement:
    def __init__(self, dict_fuzzy_crisprs, parameters, flag_dev_mode):
        self.dict_fuzzy_crisprs = dict_fuzzy_crisprs
        self.parameters = parameters
        self.flag_dev_mode = flag_dev_mode
        self.dict_fuzzy_crisprs_refined_st_end = {}
        self.dict_fuzzy_crisprs_fully_refined = {}

        self._filter_out_same_start_end_cases()
        self._filter_out_non_crispr_cases()
        self._reformat_ac_crispr_candidates()

    def _filter_out_same_start_end_cases(self):
        ssef = SameStartEndFilter(self.dict_fuzzy_crisprs)
        self.dict_fuzzy_crisprs_refined_st_end = ssef.output()

    def _filter_out_non_crispr_cases(self):
        self.param_min_avg_repeat_length = self.parameters["param_min_avg_repeat_length"]
        self.param_max_avg_repeat_length = self.parameters["param_max_avg_repeat_length"]
        self.param_min_avg_spacer_length = self.parameters["param_min_avg_spacer_length"]
        self.param_max_avg_spacer_length = self.parameters["param_max_avg_spacer_length"]
        self.param_min_repeats = self.parameters["param_min_repeats"]
        self.param_max_identical_spacers = self.parameters["param_max_identical_spacers"]
        self.param_max_identical_cluster_spacers = self.parameters["param_max_identical_cluster_spacers"]

        afsf = AdvancedFuzzySearchFilter(min_column_dominance_repeat=0.6,
                                         max_spacer_length=140, max_column_dominance_spacer=0.8,
                                         max_allowed_consecutive_spacers=self.param_max_identical_cluster_spacers,
                                         max_allowed_same_spacers=self.param_max_identical_spacers,
                                         max_inconsistent_columns=5,
                                         min_avg_repeat_length=self.param_min_avg_repeat_length,
                                         max_avg_repeat_length=self.param_max_avg_repeat_length,
                                         min_avg_spacer_length=self.param_min_avg_spacer_length,
                                         max_avg_spacer_length=self.param_max_avg_spacer_length,
                                         min_repeats=self.param_min_repeats)

        for key, values in self.dict_fuzzy_crisprs_refined_st_end.items():
            list_filtered_advanced = [afsf(value) for value in values]
            list_filtered_advanced = [x for x in list_filtered_advanced if x]
            if not list_filtered_advanced:
                sorted_by_num_errors = sorted(list(values), key=lambda x: x.number_errors)
                if sorted_by_num_errors:
                    candidate_fewer_mismatches = sorted_by_num_errors[0]
                    self.dict_fuzzy_crisprs_fully_refined[key] = [candidate_fewer_mismatches]
            else:
                self.dict_fuzzy_crisprs_fully_refined[key] = list_filtered_advanced

    def _reformat_ac_crispr_candidates(self):
        self.dict_crispr_candidates = {}
        for key, list_fuzzy in self.dict_fuzzy_crisprs_fully_refined.items():
            new_key = (key.start, key.end)
            list_crispr_candidates = [CrisprCandidate(fuzzy.list_repeats, fuzzy.list_gaped_repeats,
                                                      fuzzy.list_spacers, fuzzy.list_absolute_start)
                                      for fuzzy in list_fuzzy]

            self.dict_crispr_candidates[new_key] = list_crispr_candidates

    def output(self):
        return self.dict_crispr_candidates
