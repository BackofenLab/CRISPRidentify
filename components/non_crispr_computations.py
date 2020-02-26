import os
import subprocess
import math

from components_is_element import FullISElementSearch
from components_strand import StrandComputation
from components_cas_genes import CasGeneSearch
from components_degenerated import FullDegeneratedSearch
from components_crispr_is_merge import CrisprMerging
from components_leader_seq_search import FullLeaderSeqSearch
from components_rev_com_crispr import RevComComputation


class AdditionalCalculations:
    def __init__(self, input_file, full_dna, list_of_crisprs, hmm_model_is_element, hmm_model_cas_genes):
        self.input_file = input_file
        self.full_dna = full_dna
        self.list_of_crisprs = list_of_crisprs
        self.hmm_model_is_elements = hmm_model_is_element
        self.hmm_model_cas_genes = hmm_model_cas_genes

    def calculate_is_elements(self):
        fies = FullISElementSearch(full_dna=self.full_dna, list_of_crisprs=self.list_of_crisprs,
                                   hmm_model=self.hmm_model_is_elements, min_similarity=0.9, min_coverage=0.9)

        return fies.output()

    def calculate_cas_proteins(self):
        def _get_crispr_intervals():
            intervals = [(x.compute_stats()["start"], x.compute_stats()["end"]) for x in self.list_of_crisprs]
            return intervals

        def _filter_cas_genes(intervals, dict_cas_genes):
            dict_filtered_cas_intervals = {}
            for key, value in dict_cas_genes.items():
                for interval in intervals:
                    if interval[0] <= key[0] < interval[1]:
                        break
                    if interval[0] <= key[1] < interval[1]:
                        break
                else:
                    dict_filtered_cas_intervals[key] = value

            return dict_filtered_cas_intervals

        def _cluster_cas_genes(dict_cas_genes):
            list_clusters = []
            cluster = []
            for key in sorted(dict_cas_genes.keys()):
                value = dict_cas_genes[key]
                new_candidate = key[0], key[1], value[0], value[1]
                if not cluster:
                    cluster.append(new_candidate)
                elif abs(cluster[-1][1] - new_candidate[0]) < 500:
                    cluster.append(new_candidate)
                else:
                    list_clusters.append(cluster)
                    cluster = [new_candidate]

            if cluster:
                list_clusters.append(cluster)

            return list_clusters

        def _clusters_to_simple_representation(list_clusters):
            list_simple_clusters = []
            for cluster in list_clusters:
                cluster_start = cluster[0][0]
                cluster_end = cluster[-1][1]
                list_cas_gene_descriptions = [x[2] for x in cluster]
                list_simple_clusters.append((cluster_start, cluster_end, list_cas_gene_descriptions))
            return list_simple_clusters

        def _compute_allowed_intervals(crispr_intervals):
            allowed_interwals = []
            if not crispr_intervals:
                return [(0, math.inf)]
            else:
                allowed_interwals.append((0, crispr_intervals[0][0]))
                for index in range(len(crispr_intervals) - 1):
                    allowed_interwals.append((crispr_intervals[index][1], crispr_intervals[index+1][0]))
            allowed_interwals.append((crispr_intervals[-1][1], math.inf))
            return allowed_interwals

        def _group_by_output(allowed_intervals, list_simple_clusters):
            dict_cas_gene_order = {}
            for cluster in list_simple_clusters:
                for index, allowed_interval in enumerate(allowed_intervals):
                    if allowed_interval[0] <= cluster[0] < allowed_interval[1]:
                        if index in dict_cas_gene_order:
                            dict_cas_gene_order[index].append(cluster)
                        else:
                            dict_cas_gene_order[index] = [cluster]
                        break
            return dict_cas_gene_order

        def _group_by_output_separated(allowed_intervals, regular_clusters):
            dict_cas_gene_order_for_separated = {}
            for cluster in regular_clusters:
                for index, allowed_interval in enumerate(allowed_intervals):
                    if allowed_interval[0] <= cluster[0][0] < allowed_interval[1]:
                        if index in dict_cas_gene_order_for_separated:
                            dict_cas_gene_order_for_separated[index].append(cluster)
                        else:
                            dict_cas_gene_order_for_separated[index] = [cluster]
                        break
            return dict_cas_gene_order_for_separated

        cgs = CasGeneSearch(self.input_file, self.hmm_model_cas_genes, 0.0001)
        dict_cas_genes = cgs.dict_interval_cas
        intervals = _get_crispr_intervals()
        allowed_intervals = _compute_allowed_intervals(intervals)
        dict_filtered_cas_genes = _filter_cas_genes(intervals, dict_cas_genes)
        clustered_cas_genes = _cluster_cas_genes(dict_filtered_cas_genes)

        simple_clusters = _clusters_to_simple_representation(clustered_cas_genes)
        dict_groups = _group_by_output(allowed_intervals, simple_clusters)
        dict_groups_separated = _group_by_output_separated(allowed_intervals, clustered_cas_genes)

        return dict_groups, dict_groups_separated

    def calculate_strand(self):
        #st = StrandComputation(list_of_crisprs=self.list_of_crisprs)
        #return st.output()
        return "Forward"

    def calculate_leader(self):
        flss = FullLeaderSeqSearch(self.list_of_crisprs, self.full_dna)
        return flss.output()

    def calculate_leader_rev_com(self):
        list_of_crisprs_rev_com = [RevComComputation(crispr_candidate).output()
                                   for crispr_candidate in self.list_of_crisprs]
        flss = FullLeaderSeqSearch(list_of_crisprs_rev_com, self.full_dna)
        return flss.output()

    def calculate_degenerated(self):
        fds = FullDegeneratedSearch(full_dna=self.full_dna, list_crisprs=self.list_of_crisprs)
        return fds.output()

    def calculate_merged(self):
        fies = FullISElementSearch(full_dna=self.full_dna, list_of_crisprs=self.list_of_crisprs,
                                   hmm_model=self.hmm_model_is_elements, min_similarity=0.9, min_coverage=0.9)
        dict_is_element = fies.output()

        cm = CrisprMerging(full_dna=self.full_dna, list_crisprs=self.list_of_crisprs, dict_is_elements=dict_is_element)
        return cm.output()





