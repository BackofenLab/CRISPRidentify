import multiprocessing
from multiprocessing import Pool

from components.components_detection import VmatchRun
from components.components_detection import ClusterMaker
from components.components_detection import FilterApproximationClusters
from components.components_detection import StartEndEnhancementClusters
from components.components_detection import IntermediateEnhancementClusters
from components.components_detection import ClusterSequence
from components.components_detection import FuzzySearch


class Detection:
    def __init__(self, file_path, flags, parameters, flag_dev_mode):
        self.file_path = file_path
        self.flags = flags
        self.parameters = parameters
        self.flag_parallel = flags["flag_parallel"]
        self.flag_cpu = flags["flag_cpu"]
        self.flag_fast_run = flags["flag_fast_run"]
        self.flag_enhancement_min_max = flags["flag_enhancement_min_max"]
        self.flag_enhancement_start_end = flags["flag_enhancement_start_end"]
        self.parameters = parameters
        self.flag_dev_mode = flag_dev_mode

        self.clusters = []
        self.cluster_sequences = []
        self.dict_fuzzy_crisprs = {}

        self._get_complete_dna()
        self._run_cluster_detection()
        self._extract_cluster_sequences()
        self._run_array_detection()

    def _get_complete_dna(self):
        with open(self.file_path, 'r') as f:
            lines = f.readlines()

        self.input_header = lines[0]
        self.dna = ''.join([line.strip() for line in lines[1:]])
        self.dna_length = len(self.dna)
        self.dna = self.dna.upper()

    def _run_cluster_detection(self):
        vr = VmatchRun(self.file_path, self.flag_fast_run)
        list_repeats_from_vmatch = vr.output()
        #print("list vmatch repeats", list_repeats_from_vmatch)

        cm = ClusterMaker(list_repeats_from_vmatch, self.dna)
        self.clusters = cm.output()

        fa = FilterApproximationClusters(self.clusters)
        self.clusters = fa.output()

        st = StartEndEnhancementClusters(self.clusters)
        self.clusters = st.output()

        ie = IntermediateEnhancementClusters(self.clusters)
        self.clusters = ie.output()

    def _extract_cluster_sequences(self):
        for cluster in self.clusters:
            seq_start = max(0, cluster.begin - 100)
            seq_end = min(len(self.dna), cluster.end + 100)
            cluster_seq = self.dna[seq_start:seq_end]
            tup_cluster_dif_rep = tuple(cluster.list_clust_dif_rep_seq)

            self.cluster_sequences.append(ClusterSequence(cluster_seq, seq_start, seq_end, tup_cluster_dif_rep))

    @staticmethod
    def _parallel_run_fuzzy_run(input_tuple):
        repeat, sequence, start, weighted_error = input_tuple

        return FuzzySearch(sequence, start,
                           repeat, weighted_error)

    def _run_array_detection(self):
        weighted_error = "{i<=3,d<=3,s<=3,i+d+s<=6}"
        parallel = self.flag_parallel

        if parallel:
            for cluster_sequence in self.cluster_sequences:
                nr = len(cluster_sequence.tuple_repeats)
                input_tuples = zip(cluster_sequence.tuple_repeats, [cluster_sequence.sequence] * nr,
                                   [cluster_sequence.start] * nr, [weighted_error] * nr)

                num_workers_suggested = multiprocessing.cpu_count() if self.flag_cpu == "ALL" else int(self.flag_cpu)
                max_possible = multiprocessing.cpu_count()
                num_workers = num_workers_suggested if num_workers_suggested < max_possible else max_possible
                with Pool(num_workers) as p:
                    fuzzy_results = p.map(self._parallel_run_fuzzy_run, input_tuples)
                    fuzzy_results = [x for x in fuzzy_results if x.match_hit]
                    fuzzy_results = [x for x in fuzzy_results if len(x.list_repeats) > 1]

                self.dict_fuzzy_crisprs[cluster_sequence] = fuzzy_results
        else:
            for cluster_sequence in self.cluster_sequences:
                list_fuzzy_results = []
                for repeat in cluster_sequence.tuple_repeats:
                    fuzzy_s = FuzzySearch(cluster_sequence.sequence, cluster_sequence.start,
                                          repeat, weighted_error)
                    if fuzzy_s.match_hit:
                        if len(fuzzy_s.list_repeats) > 1:
                            list_fuzzy_results.append(fuzzy_s)

                self.dict_fuzzy_crisprs[cluster_sequence] = list_fuzzy_results

    def output(self):
        return self.dict_fuzzy_crisprs