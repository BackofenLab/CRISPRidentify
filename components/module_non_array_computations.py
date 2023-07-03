import math

from components.components_non_array_computations import StrandComputation
from components.components_non_array_computations import StrandComputationNew
from components.components_non_array_computations import FullISElementSearch
from components.components_non_array_computations import complete_info_with_cas_identifier
from components.components_non_array_computations import FullLeaderSeqSearch
from components.components_non_array_computations import RevComComputation


class NonArrayComputations:
    def __init__(self, file_path, categories, flags_non_arrays_computations, flag_dev_mode, absolute_directory_path):
        self.file_path = file_path
        self.categories = categories
        self.flags_non_arrays_computations = flags_non_arrays_computations
        self.flag_dev_mode=flag_dev_mode
        self.absolute_directory_path = absolute_directory_path

        self.list_of_crisprs_bona_fide = [self.categories[0][key][0][1] for key in sorted(self.categories[0].keys())]
        self.list_of_crisprs_alternative = [el[1] for key in self.categories[1].keys()
                                            for el in self.categories[1][key]]
        self.list_of_crisprs_possible = [el[1] for key in self.categories[2].keys()
                                         for el in self.categories[2][key]]

        self.hmm_model_is_elements = "tools/hmm_search/models_is_element.hmm"

        self.is_element_result = {}
        self.cas_results = {}
        self.cassete_results = {}
        self.unstructured_cas_result_from_cas_identifier = {}
        self.strand_results = {}
        self.leader_results = {}
        self.downstream_results = {}
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

        if self.flags_non_arrays_computations["flag_cas"]:
            self._calculate_cas_proteins()
        if self.flags_non_arrays_computations["flag_is"]:
            self._calculate_is_elements()

        self.data_with_all_computations = {"IS": self.is_element_result,
                                           "Cas": self.cas_results,
                                           "Strand": self.strand_results,
                                           "Leader": [self.leader_results_bona_fide, self.leader_results_alternative, self.leader_results_possible],
                                           "Downstream": [self.downstream_results_bona_fide, self.downstream_results_alternative, self.downstream_results_possible],
                                           "Unstructured_Cas":self.unstructured_cas_result_from_cas_identifier,
                                           "Cassettes": self.cassete_results}

    def _calculate_is_elements(self):
        fies = FullISElementSearch(full_dna=self.dna, list_of_crisprs=self.list_of_crisprs_bona_fide,
                                   hmm_model=self.hmm_model_is_elements, min_similarity=0.9, min_coverage=0.9)

        self.is_element_result = fies.output()

    def _calculate_cas_proteins(self):
        def _get_crispr_intervals():
            intervals = [(x.compute_stats()["start"], x.compute_stats()["end"]) for x in self.list_of_crisprs_bona_fide]
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

        dict_cas_genes, dict_cassete_labels = complete_info_with_cas_identifier(self.file_path,
                                                                                self.absolute_directory_path)

        self.cassete_results = dict_cassete_labels
        self.unstructured_cas_result_from_cas_identifier = dict_cas_genes

        intervals = _get_crispr_intervals()
        allowed_intervals = _compute_allowed_intervals(intervals)
        dict_filtered_cas_genes = _filter_cas_genes(intervals, dict_cas_genes)
        clustered_cas_genes = _cluster_cas_genes(dict_filtered_cas_genes)

        simple_clusters = _clusters_to_simple_representation(clustered_cas_genes)
        dict_groups = _group_by_output(allowed_intervals, simple_clusters)
        #dict_groups_separated = _group_by_output_separated(allowed_intervals, clustered_cas_genes)

        self.cas_results = dict_groups

    def _calculate_strand(self):
        if self.flags_non_arrays_computations["flag_strand"]:
            st = StrandComputationNew(list_of_crisprs=self.list_of_crisprs_bona_fide,
                                      absolute_directory_path=self.absolute_directory_path)
            self.strand_results["Bona-fide"] = st.output()
            st = StrandComputationNew(list_of_crisprs=self.list_of_crisprs_alternative,
                                      absolute_directory_path=self.absolute_directory_path)
            self.strand_results["Alternative"] = st.output()
            st = StrandComputationNew(list_of_crisprs=self.list_of_crisprs_possible,
                                      absolute_directory_path=self.absolute_directory_path)
            self.strand_results["Possible"] = st.output()


            #except Exception:
            #    st = StrandComputation(list_of_crisprs=self.list_of_crisprs_bona_fide,
            #                           absolute_directory_path=self.absolute_directory_path)
            #    self.strand_results["Bona-fide"] = st.output()
            #    st = StrandComputation(list_of_crisprs=self.list_of_crisprs_alternative,
            #                           absolute_directory_path=self.absolute_directory_path)
            #    self.strand_results["Alternative"] = st.output()
            #    st = StrandComputation(list_of_crisprs=self.list_of_crisprs_possible,
            #                           absolute_directory_path=self.absolute_directory_path)
            #    self.strand_results["Possible"] = st.output()
        else:
            self.strand_results["Bona-fide"] = {index: "Forward (Orientation was not computed)"
                                                for index in range(len(self.list_of_crisprs_bona_fide))}
            self.strand_results["Alternative"] = {index: "Forward (Orientation was not computed)"
                                                  for index in range(len(self.list_of_crisprs_alternative))}
            self.strand_results["Possible"] = {index: "Forward (Orientation was not computed)"
                                               for index in range(len(self.list_of_crisprs_possible))}

    def _calculate_leader(self):
        flss_bona_fide = FullLeaderSeqSearch(self.list_of_crisprs_bona_fide, self.strand_results["Bona-fide"], self.dna)
        self.leader_results_bona_fide, self.downstream_results_bona_fide = flss_bona_fide.output()

        flss_alternative = FullLeaderSeqSearch(self.list_of_crisprs_alternative, self.strand_results["Alternative"],
                                               self.dna)
        self.leader_results_alternative, self.downstream_results_alternative = flss_alternative.output()

        flss_possible = FullLeaderSeqSearch(self.list_of_crisprs_possible, self.strand_results["Possible"], self.dna)
        self.leader_results_possible, self.downstream_results_possible = flss_possible.output()

    def output(self):
        return self.data_with_all_computations
