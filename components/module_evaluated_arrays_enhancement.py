from components_evaluated_arrays_enhancement import IterativeDegeneratedSearch
from components_evaluated_arrays_enhancement import create_boundaries_for_intervals
from components_evaluated_arrays_enhancement import ArrayRefinerInsertionsDeletions


class EvaluatedArraysEnhancement:
    def __init__(self, file_path, categories, parameters):
        self.file_path = file_path
        self.categories = categories
        self.parameters = parameters

        self.bona_fide_arrays = categories[0]
        self.possible_arrays = categories[2]

        self.dict_arrays_into_categories_enhanced = {}

        self._get_complete_dna()
        self._search_missed_or_degenerated_repeats()
        self._refine_nucleotides_repeat_spacer()

    def _get_complete_dna(self):
        with open(self.file_path, 'r') as f:
            lines = f.readlines()

        self.input_header = lines[0]
        self.dna = ''.join([line.strip() for line in lines[1:]])
        self.dna_length = len(self.dna)
        self.dna = self.dna.upper()

    def _search_missed_or_degenerated_repeats(self):
        intervals = []
        arrays = []
        for interval, list_data in self.bona_fide_arrays.items():
            intervals.append(interval)
            arrays.append(list_data[0][1])

        boundaries = create_boundaries_for_intervals(intervals, 500)

        for interval, array, boundary in zip(intervals, arrays, boundaries):

            consensus = array.consensus
            list_repeats = array.list_repeats
            list_repeats_starts = array.list_repeat_starts
            list_spacers = array.list_spacers

            ids = IterativeDegeneratedSearch(full_dna=self.dna,
                                             repeat_seq_candidate=consensus,
                                             spacer_margin=self.parameters["param_spacer_margin_degenerated_search"],
                                             repeat_seq_candidate_gaped=None,
                                             list_repeats_starts=list_repeats_starts,
                                             list_repeats=list_repeats,
                                             list_spacers=list_spacers,
                                             start_flanking_region_left=boundary[0],
                                             end_flanking_region_right=boundary[1],
                                             allowed_max_editing_distance=self.parameters["param_max_edit_distance"],
                                             iterative_size_flanking_region=150,
                                             prevent_long_spacers=True,
                                             attempt_to_improve_initial_array=True)

            new_crispr_candidate = ids.output()
            self.bona_fide_arrays[interval][0][1] = new_crispr_candidate

    def _refine_nucleotides_repeat_spacer(self):
        intervals = []
        arrays = []
        for interval, list_data in self.bona_fide_arrays.items():
            intervals.append(interval)
            arrays.append(list_data[0][1])

        for interval, array in zip(intervals, arrays):
            arid = ArrayRefinerInsertionsDeletions(array)
            new_crispr_candidate = arid.output()
            self.bona_fide_arrays[interval][0][1] = new_crispr_candidate

    def output(self):
        self.categories[0] = self.bona_fide_arrays
        return self.categories