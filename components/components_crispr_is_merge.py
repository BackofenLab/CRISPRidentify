import os
from feature_extraction import CrisprSimilarity
from detection import FuzzySearch
from crispr_candidate import CrisprCandidate


class CrisprMerging:
    def __init__(self, full_dna, list_crisprs, dict_is_elements):
        self.full_dna = full_dna
        self.list_crisprs = list_crisprs
        self.dict_is_elements = dict_is_elements

        self.clusters = []

        self._cluster_search()
        self._new_crisprs_creation()

    def _cluster_search(self):
        dict_is = self.dict_is_elements
        list_index_is = list(dict_is.keys())

        clusters = []
        current_cluster = []
        for index, crispr_array in enumerate(self.list_crisprs):
            if not current_cluster:
                current_cluster.append(index)
            else:
                prev_index = index - 1
                prev_crispr = self.list_crisprs[prev_index]
                if prev_index in list_index_is:
                    if CheckMerging(prev_crispr, crispr_array, self.full_dna, 3000).check_merging():
                        current_cluster.append(index)
                    else:
                        clusters.append(current_cluster)
                        current_cluster = [index]
                else:
                    if CheckMerging(prev_crispr, crispr_array, self.full_dna, 500).check_merging():
                        current_cluster.append(index)
                    else:
                        clusters.append(current_cluster)
                        current_cluster = [index]

        clusters.append(current_cluster)
        self.clusters = [cluster for cluster in clusters if len(cluster) > 1]

    def _new_crisprs_creation(self):
        self.new_crisprs = []
        for cluster in self.clusters:
            first_crispr_index = cluster[0]
            last_crispr_index = cluster[-1]
            first_crispr = self.list_crisprs[first_crispr_index]
            last_crispr = self.list_crisprs[last_crispr_index]

            cluster_start = first_crispr.list_repeat_starts[0]
            cluster_end = last_crispr.list_repeat_starts[-1] + len(last_crispr.list_repeats[-1])
            consensus = first_crispr.consensus
            interval = self.full_dna[cluster_start:cluster_end]
            fuzzy_search = FuzzySearch(sequence=interval, sequence_start=cluster_start,
                                       repeat=consensus, weighted_error="{i<=3,d<=3,s<=3,i+d+s<=6}")
            cc = CrisprCandidate(fuzzy_search.list_repeats, fuzzy_search.list_gaped_repeats,
                                 fuzzy_search.list_spacers, fuzzy_search.list_absolute_start)
            self.new_crisprs.append(cc)

    def output(self):
        return {index: cc for index, cc in enumerate(self.new_crisprs)}


class CheckMerging:
    def __init__(self, crispr_1, crispr_2, full_dna, max_distance):
        self.crispr_1 = crispr_1
        self.crispr_2 = crispr_2
        self.full_dna = full_dna
        self.max_distance = max_distance

    def check_merging(self):
        if self._check_merging_consensus():
            if self._check_merging_intervals():
                if self._check_merging_flanks():
                    return True
        return False

    def _check_merging_consensus(self):
        consensus_1 = self.crispr_1.consensus
        consensus_2 = self.crispr_2.consensus
        if consensus_1 == consensus_2:
            return True
        return False

    def _check_merging_intervals(self):
        self.crispr_1_start = self.crispr_1.list_repeat_starts[0]
        self.crispr_1_end = self.crispr_1.list_repeat_starts[-1] + len(self.crispr_1.list_repeats[-1])

        self.crispr_2_start = self.crispr_2.list_repeat_starts[0]
        self.crispr_2_end = self.crispr_2.list_repeat_starts[-1] + len(self.crispr_2.list_repeats[-1])
        if abs(self.crispr_1_end - self.crispr_2_start) <= self.max_distance:
            return True
        return False

    def _check_merging_flanks(self):
        flank_1 = self.full_dna[self.crispr_1_start - 100:self.crispr_1_start]
        flank_2 = self.full_dna[self.crispr_1_end:self.crispr_1_end+100]

        flank_3 = self.full_dna[self.crispr_2_start - 100:self.crispr_2_start]
        flank_4 = self.full_dna[self.crispr_2_end:self.crispr_2_end + 100]

        combination_1 = CrisprSimilarity(0, [flank_1, flank_3], []).output()[0]
        combination_2 = CrisprSimilarity(0, [flank_1, flank_4], []).output()[0]

        combination_3 = CrisprSimilarity(0, [flank_2, flank_3], []).output()[0]
        combination_4 = CrisprSimilarity(0, [flank_2, flank_4], []).output()[0]

        if max(combination_1, combination_2, combination_3, combination_4) > 0.7:
            return False
        return True








