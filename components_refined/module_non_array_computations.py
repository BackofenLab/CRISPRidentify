import math

from components_non_array_computations import StrandComputation
from components_non_array_computations import FullISElementSearch
from components_non_array_computations import complete_info_with_cas_identifier
from components_non_array_computations import FullLeaderSeqSearch
from components_non_array_computations import RevComComputation


class NonArrayComputations:
    def __init__(self, file_path, categories, flags_non_arrays_computations):
        self.file_path = file_path
        self.categories = categories
        self.flags_non_arrays_computations = flags_non_arrays_computations

        self.list_of_crisprs = [list_data[0][1] for list_data in self.categories[0].values()]
        self.hmm_model_is_elements = "tools/hmm_search/models_is_element.hmm"

        self.is_element_result = {}
        self.cas_results = {}
        self.strand_results = {}
        self.leader_results = {}
        self.leader_results_rev_com = {}
        self.data_with_all_computations = {}

        self._get_complete_dna()
        self._calculate_all_non_array_values()

    def _get_complete_dna(self):
        with open(self.file_path, 'r') as f:
            lines = f.readlines()

        self.input_header = lines[0]
        self.dna = ''.join([line.strip() for line in lines[1:]])
        self.dna_length = len(self.dna)
        self.dna = self.dna.upper()

    def _calculate_all_non_array_values(self):
        self._calculate_strand()
        self._calculate_leader()
        self._calculate_leader_rev_com()

        if self.flags_non_arrays_computations["flag_cas"]:
            self._calculate_cas_proteins()
        if self.flags_non_arrays_computations["flag_is"]:
            self._calculate_is_elements()

        self.data_with_all_computations = {"IS": self.is_element_result,
                                           "Cas": self.cas_results,
                                           "Strand": self.strand_results,
                                           "Leader": self.leader_results,
                                           "Leader_rev_com": self.leader_results_rev_com}

    def _calculate_is_elements(self):
        fies = FullISElementSearch(full_dna=self.dna, list_of_crisprs=self.list_of_crisprs,
                                   hmm_model=self.hmm_model_is_elements, min_similarity=0.9, min_coverage=0.9)

        self.is_element_result = fies.output()

    def _calculate_cas_proteins(self):
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
                new_candidate = key[0], key[1], value
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

        dict_cas_genes = complete_info_with_cas_identifier(self.file_path)

        intervals = _get_crispr_intervals()
        allowed_intervals = _compute_allowed_intervals(intervals)
        dict_filtered_cas_genes = _filter_cas_genes(intervals, dict_cas_genes)
        clustered_cas_genes = _cluster_cas_genes(dict_filtered_cas_genes)

        simple_clusters = _clusters_to_simple_representation(clustered_cas_genes)
        dict_groups = _group_by_output(allowed_intervals, simple_clusters)
        #dict_groups_separated = _group_by_output_separated(allowed_intervals, clustered_cas_genes)

        self.cas_results = dict_groups

    def _calculate_strand(self):
        st = StrandComputation(list_of_crisprs=self.list_of_crisprs)
        self.strand_results = st.output()

    def _calculate_leader(self):
        flss = FullLeaderSeqSearch(self.list_of_crisprs, self.dna)
        self.leader_results = flss.output()

    def _calculate_leader_rev_com(self):
        list_of_crisprs_rev_com = [RevComComputation(crispr_candidate).output()
                                   for crispr_candidate in self.list_of_crisprs]
        flss = FullLeaderSeqSearch(list_of_crisprs_rev_com, self.dna)
        self.leader_results_rev_com = flss.output()

    def output(self):
        return self.data_with_all_computations