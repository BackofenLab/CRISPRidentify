import regex

from components.components_detection_refinement import CrisprCandidate


class IterativeDegeneratedSearch:
    def __init__(self, full_dna, repeat_seq_candidate,
                 spacer_margin, repeat_seq_candidate_gaped,
                 list_repeats_starts, list_repeats, list_spacers,
                 start_flanking_region_left, end_flanking_region_right,
                 allowed_max_editing_distance,
                 iterative_size_flanking_region, prevent_long_spacers=True,
                 attempt_to_improve_initial_array=True):

        self.final_coordinates = None
        self.full_dna = full_dna
        self.repeat_seq_candidate = repeat_seq_candidate
        self.spacer_margin = spacer_margin
        self.repeats_seq_candidate_gaped = repeat_seq_candidate_gaped
        self.list_repeats_starts = list_repeats_starts
        self.list_repeats = list_repeats
        self.list_spacers = list_spacers
        self.start_flanking_region_left = start_flanking_region_left
        self.end_flanking_region_right = end_flanking_region_right
        self.iterative_size_flanking_region = iterative_size_flanking_region
        self.allowed_max_editing_distance = allowed_max_editing_distance
        self.prevent_long_spacers = prevent_long_spacers
        self.attempt_to_improve_initial_array = attempt_to_improve_initial_array

        self.flag_make_step_left = True
        self.flag_make_step_right = True

        self._calculate_representation_original_array()
        self._calculate_max_possible_new_spacer_length()

        self._left_flank_iterative_search()
        self._right_flank_iterative_search()
        self._build_new_representation()

    def _calculate_representation_original_array(self):
        original_interval_start = max(0, self.list_repeats_starts[0] - 5)  # flanks for regex
        original_interval_end = self.list_repeats_starts[-1] + len(self.list_repeats[-1]) + 5  # flanks for regex

        original_dna_interval = self.full_dna[original_interval_start:original_interval_end]
        search_pattern = f"(?e)({self.repeat_seq_candidate})" + "{" + "i<=3,d<=3,s<=3,s+i+d<=6" + "}"

        original_matches = list(regex.finditer(search_pattern, original_dna_interval))

        self.list_repeats_original = [match.group() for match in original_matches]

        self.list_repeats_starts_original = [original_interval_start + match.start() for match in original_matches]

        if len(original_matches) > 1:
            internal_spacer_starts = [match.end() for match in original_matches][:-1]
            internal_spacer_end = [match.start() for match in original_matches][1:]
            internal_spacer_coordinates = [(start, end) for start, end in
                                           zip(internal_spacer_starts, internal_spacer_end)]
            self.list_spacers_original = [original_dna_interval[start_end[0]:start_end[1]]
                                          for start_end in internal_spacer_coordinates]
        else:
            self.list_spacers_original = []

        relative_error_indexes = []
        for match in original_matches:
            tuple_match_errors = match.fuzzy_changes
            list_relative_errors = [[e - match.start() for e in err_type] for err_type in tuple_match_errors]
            relative_error_indexes.append(list_relative_errors)

        self.list_repeats_error_indexes_original = relative_error_indexes

        self.list_fuzzy_counts_original = [match.fuzzy_counts for match in original_matches]

    def _calculate_max_possible_new_spacer_length(self):
        if self.list_spacers:
            self.max_possible_new_spacer_len = max((len(x) for x in self.list_spacers)) + self.spacer_margin
        else:
            self.max_possible_new_spacer_len = 60

    def _make_new_step_left(self, interval_start):
        if interval_start is not None:
            interval_end = interval_start
            interval_start = max(self.start_flanking_region_left,
                                 interval_start - self.iterative_size_flanking_region)

            interval_dna = self.full_dna[interval_start:interval_end]
            distance = 6
            return interval_start, interval_end, interval_dna, distance
        else:
            interval_start = max(self.start_flanking_region_left,
                                 self.list_repeats_starts[0] - self.iterative_size_flanking_region)
            interval_end = self.list_repeats_starts[0]

            interval_dna = self.full_dna[interval_start:interval_end]
            distance = 6
            return interval_start, interval_end, interval_dna, distance

    def _make_new_step_right(self, interval_end):
        if interval_end is not None:
            interval_start = interval_end
            interval_end = min(self.end_flanking_region_right,
                               interval_end + self.iterative_size_flanking_region)

            interval_dna = self.full_dna[interval_start:interval_end]
            distance = 6
            return interval_start, interval_end, interval_dna, distance
        else:
            interval_start = self.list_repeats_starts[-1] + len(self.list_repeats[-1])
            interval_end = min(self.end_flanking_region_right,
                               interval_start + self.iterative_size_flanking_region)

            interval_dna = self.full_dna[interval_start:interval_end]
            distance = 6
            return interval_start, interval_end, interval_dna, distance

    def increase_editing_distance(self, distance):
        if distance < self.allowed_max_editing_distance:
            return distance + 2
        else:
            return self.allowed_max_editing_distance

    def _create_pattern(self, distance):
        half_distance = int(distance / 2)
        pattern = f"(?be)({self.repeat_seq_candidate})" + "{" + \
                  f"s<={half_distance},i<={half_distance},d<={half_distance},s+i+d<={distance}" + "}"
        return pattern

    def _left_flank_iterative_search(self):
        self.list_repeats_in_intervals_left = []
        self.list_repeats_starts_in_intervals_left = []
        self.list_repeats_errors_in_intervals_left = []
        self.list_fuzzy_counts_in_intervals_left = []
        self.list_spacers_in_intervals_left = []
        self.list_bridge_spacers_left = []

        flag_stop_iterations = False
        most_left_repeat_match_start = None

        while not flag_stop_iterations:
            if self.flag_make_step_left:
                new_step_left = self._make_new_step_left(most_left_repeat_match_start)
                interval_start, interval_end, interval_dna, distance = new_step_left
            else:
                distance = self.increase_editing_distance(distance)

            search_pattern = self._create_pattern(distance)
            ibfs = IntervalBasedFuzzySearch(interval_start, interval_end, interval_dna, search_pattern)

            if ibfs.flag_exist_match:
                bridge_spacer = ibfs.compute_left_bridge_spacer()
                len_bridge_spacer_to_the_left = len(bridge_spacer)
                if len_bridge_spacer_to_the_left < self.max_possible_new_spacer_len:
                    list_repeats = ibfs.get_repeats()
                    list_spacers = ibfs.get_internal_spacers()
                    list_repeats_errors = ibfs.get_relative_error_indexes()

                    list_repeats_relative_starts = ibfs.get_relative_starts()
                    list_repeats_starts = [(x + interval_start) for x in list_repeats_relative_starts]
                    most_left_repeat_match_start = list_repeats_starts[0]
                    list_fuzzy_counts = ibfs.get_errors()

                    self.list_repeats_in_intervals_left.append(list_repeats)
                    self.list_repeats_starts_in_intervals_left.append(list_repeats_starts)
                    self.list_repeats_errors_in_intervals_left.append(list_repeats_errors)
                    self.list_fuzzy_counts_in_intervals_left.append(list_fuzzy_counts)
                    self.list_spacers_in_intervals_left.append(list_spacers)
                    self.list_bridge_spacers_left.append(bridge_spacer)

                    self.flag_make_step_left = True
                else:
                    self.flag_make_step_left = False
                    if distance == self.allowed_max_editing_distance:
                        flag_stop_iterations = True
            else:
                self.flag_make_step_left = False
                if distance == self.allowed_max_editing_distance:
                    flag_stop_iterations = True

    def _right_flank_iterative_search(self):
        self.list_repeats_in_intervals_right = []
        self.list_repeats_starts_in_intervals_right = []
        self.list_repeats_errors_in_intervals_right = []
        self.list_fuzzy_counts_in_intervals_right = []
        self.list_spacers_in_intervals_right = []
        self.list_bridge_spacers_right = []

        flag_stop_iterations = False
        most_right_repeat_match_end = None

        while not flag_stop_iterations:
            if self.flag_make_step_right:
                new_step_right = self._make_new_step_right(most_right_repeat_match_end)
                interval_start, interval_end, interval_dna, distance = new_step_right
            else:
                distance = self.increase_editing_distance(distance)

            search_pattern = self._create_pattern(distance)
            ibfs = IntervalBasedFuzzySearch(interval_start, interval_end, interval_dna, search_pattern)

            if ibfs.flag_exist_match:
                bridge_spacer = ibfs.compute_right_bridge_spacer()
                len_bridge_spacer_to_the_right = len(bridge_spacer)
                if len_bridge_spacer_to_the_right < self.max_possible_new_spacer_len:
                    list_repeats = ibfs.get_repeats()
                    list_spacers = ibfs.get_internal_spacers()
                    list_repeats_errors = ibfs.get_relative_error_indexes()

                    list_repeats_relative_starts = ibfs.get_relative_starts()
                    list_repeats_starts = [(x + interval_start) for x in list_repeats_relative_starts]
                    most_right_repeat_match_end = list_repeats_starts[-1] + len(list_repeats[-1])
                    list_fuzzy_counts = ibfs.get_errors()

                    self.list_repeats_in_intervals_right.append(list_repeats)
                    self.list_repeats_starts_in_intervals_right.append(list_repeats_starts)
                    self.list_repeats_errors_in_intervals_right.append(list_repeats_errors)
                    self.list_fuzzy_counts_in_intervals_right.append(list_fuzzy_counts)
                    self.list_spacers_in_intervals_right.append(list_spacers)
                    self.list_bridge_spacers_right.append(bridge_spacer)

                    self.flag_make_step_right = True
                else:
                    self.flag_make_step_right = False
                    if distance == self.allowed_max_editing_distance:
                        flag_stop_iterations = True
            else:
                self.flag_make_step_right = False
                if distance == self.allowed_max_editing_distance:
                    flag_stop_iterations = True

    def _build_new_representation(self):
        self.list_repeats_in_intervals_left = self.list_repeats_in_intervals_left[::-1]
        self.list_repeats_starts_in_intervals_left = self.list_repeats_starts_in_intervals_left[::-1]
        self.list_repeats_errors_in_intervals_left = self.list_repeats_errors_in_intervals_left[::-1]
        self.list_fuzzy_counts_in_intervals_left = self.list_fuzzy_counts_in_intervals_left[::-1]
        self.list_spacers_in_intervals_left = self.list_spacers_in_intervals_left[::-1]
        self.list_bridge_spacers_left = self.list_bridge_spacers_left[::-1]

        list_repeats_left = [r for list_rep in self.list_repeats_in_intervals_left for r in list_rep]
        list_repeats_starts_left = [s for list_starts in self.list_repeats_starts_in_intervals_left for s in
                                    list_starts]
        list_repeats_errors_in_interval_left = [l_err for list_all_errors in self.list_repeats_errors_in_intervals_left
                                                for l_err in list_all_errors]

        list_fuzzy_counts_left = [fuzzy_count for list_fc in self.list_fuzzy_counts_in_intervals_left
                                  for fuzzy_count in list_fc]

        list_spacers_left = []
        for internal_spacers, bridge_spacer in zip(self.list_spacers_in_intervals_left, self.list_bridge_spacers_left):
            list_spacers_left += internal_spacers
            list_spacers_left.append(bridge_spacer)

        list_repeats_right = [r for list_rep in self.list_repeats_in_intervals_right for r in list_rep]
        list_repeats_starts_right = [s for list_starts in self.list_repeats_starts_in_intervals_right for s in
                                     list_starts]
        list_repeats_errors_in_interval_right = [l_err for list_all_errors in
                                                 self.list_repeats_errors_in_intervals_right
                                                 for l_err in list_all_errors]

        list_fuzzy_counts_right = [fuzzy_count for list_fc in self.list_fuzzy_counts_in_intervals_right
                                   for fuzzy_count in list_fc]

        list_spacers_right = []
        for internal_spacers, bridge_spacer in zip(self.list_spacers_in_intervals_right,
                                                   self.list_bridge_spacers_right):
            list_spacers_right += internal_spacers
            list_spacers_right.append(bridge_spacer)

        list_repeats = list_repeats_left + self.list_repeats_original + list_repeats_right
        list_repeats_starts = list_repeats_starts_left + self.list_repeats_starts_original + list_repeats_starts_right
        list_spacers = list_spacers_left + self.list_spacers_original + list_spacers_right

        list_repeats_error_indexes = list_repeats_errors_in_interval_left + self.list_repeats_error_indexes_original + \
                                     list_repeats_errors_in_interval_right

        list_fuzzy_counts = list_fuzzy_counts_left + self.list_fuzzy_counts_original + list_fuzzy_counts_right

        self.final_coordinates = list_repeats_starts[0], list_repeats_starts[-1] + len(list_repeats[-1])

        dot_representation_maker = DotRepresentationMaker(list_repeats=list_repeats,
                                                          list_repeats_starts=list_repeats_starts,
                                                          list_spacers=list_spacers,
                                                          consensus=self.repeat_seq_candidate,
                                                          relative_error_indexes=list_repeats_error_indexes,
                                                          fuzzy_counts=list_fuzzy_counts)

        list_repeats_gaped = dot_representation_maker.list_gaped_repeats

        self.crispr_candidate = CrisprCandidate(list_repeats, list_repeats_gaped,
                                                list_spacers, list_repeats_starts)

    def output(self):
        return self.crispr_candidate


class IntervalBasedFuzzySearch:
    def __init__(self, interval_start, interval_end, interval_dna, search_pattern):
        self.interval_start = interval_start
        self.interval_end = interval_end
        self.interval_dna = interval_dna
        self.search_pattern = search_pattern

        self.len_bridge_spacer_left = None
        self.len_bridge_spacer_right = None

        self._find_match()

    def _find_match(self):
        self.matches = list(regex.finditer(self.search_pattern, self.interval_dna))
        if self.matches:
            self.flag_exist_match = True
        else:
            self.flag_exist_match = False

    def compute_left_bridge_spacer(self):
        bridge_spacer_start = self.matches[-1].end()
        bridge_spacer = self.interval_dna[bridge_spacer_start:]
        return bridge_spacer

    def compute_right_bridge_spacer(self):
        bridge_spacer_end = self.matches[0].start()
        bridge_spacer = self.interval_dna[:bridge_spacer_end]
        return bridge_spacer

    def get_repeats(self):
        list_repeats = [match.group() for match in self.matches]
        return list_repeats

    def get_internal_spacers(self):
        if len(self.matches) > 1:
            internal_spacer_starts = [match.end() for match in self.matches][:-1]
            internal_spacer_end = [match.start() for match in self.matches][1:]
            internal_spacer_coordinates = [(start, end) for start, end in
                                           zip(internal_spacer_starts, internal_spacer_end)]
            list_internal_spacers = [self.interval_dna[start_end[0]:start_end[1]]
                                     for start_end in internal_spacer_coordinates]
        else:
            list_internal_spacers = []

        return list_internal_spacers

    def get_errors(self):
        return [match.fuzzy_counts for match in self.matches]

    def get_relative_starts(self):
        return [match.start() for match in self.matches]

    def get_relative_error_indexes(self):
        relative_error_indexes = []
        for match in self.matches:
            tuple_match_errors = match.fuzzy_changes
            list_relative_errors = [[e - match.start() for e in err_type] for err_type in tuple_match_errors]
            relative_error_indexes.append(list_relative_errors)

        return relative_error_indexes


def get_full_dna(dna_folder, acc_num):
    with open(f"{dna_folder}/{acc_num}.fasta") as f:
        lines = f.readlines()

    dna = "".join([l.strip() for l in lines[1:]])
    return dna


def get_all_candidates_intervals(pickle_result):
    best_results = pickle_result["best"]
    intervals = list(best_results.keys())
    return sorted(intervals)


class DotRepresentationMaker:
    def __init__(self, list_repeats, list_repeats_starts, list_spacers, consensus,
                 relative_error_indexes, fuzzy_counts):
        self.list_repeats = list_repeats
        self.list_spacers = list_spacers
        self.list_relative_error_indexes = relative_error_indexes
        self.consensus = consensus
        self.list_repeat_absolute_starts = list_repeats_starts
        self.fuzzy_counts = fuzzy_counts

        self.list_all_insertions = []
        self.list_all_deletions = []

        self._create_representation()
        self._get_transitional_repeat()

    def _create_representation(self):
        blocks = []
        for repeat, repeat_relative_errors in zip(self.list_repeats, self.list_relative_error_indexes):
            gaped_consensus, new_repeat_representation = search_pair_handler(repeat, self.consensus,
                                                                             repeat_relative_errors)


            block = BlockRepeats(gaped_consensus, [new_repeat_representation])
            blocks.append(block)

        br = BlockRecomputation(blocks)
        self.gaped_matching_pattern, self.list_gaped_repeats = br.output_consensus_and_repeats()
        self.list_dotted_gaped_repeats = [self.apply_dots(x, self.gaped_matching_pattern)
                                          for x in self.list_gaped_repeats]

    def _get_transitional_repeat(self):
        dict_transitions = {frozenset(("A", "T")): "W",
                            frozenset(("C", "G")): "S",
                            frozenset(("A", "C")): "M",
                            frozenset(("G", "T")): "K",
                            frozenset(("A", "G")): "R",
                            frozenset(("C", "T")): "Y"}

        transition_repeat = []
        for index, char in enumerate(self.gaped_matching_pattern):
            column = [r[index] for r in self.list_gaped_repeats]
            column = set(column)
            difference = column - {char}
            if len(difference) == 1:
                alternative_char = list(difference)[0]
                if frozenset((alternative_char, char)) in dict_transitions:
                    transition_char = dict_transitions[frozenset((alternative_char, char))]
                    transition_repeat.append(transition_char)
                else:
                    transition_repeat.append(char)
            else:
                transition_repeat.append(char)

        transition_repeat = "".join(transition_repeat)
        transition_repeat = transition_repeat.replace(" ", "")
        self.transition_consensus = transition_repeat

    @staticmethod
    def apply_dots(repeat_with_gaps, matching_pattern_with_gaps):
        dotted_repeat = ''
        if len(repeat_with_gaps) != len(matching_pattern_with_gaps):
            raise ValueError
        for char_repeat, char_pattern in zip(repeat_with_gaps, matching_pattern_with_gaps):
            if char_repeat == char_pattern:
                if char_repeat == ' ':
                    dotted_repeat += char_repeat
                else:
                    dotted_repeat += '.'
            else:
                dotted_repeat += char_repeat

        return dotted_repeat

    def create_dot_representation(self):
        dot_repr = ""
        max_len_repeat_start = max(len(str(x)) for x in self.list_repeat_absolute_starts)
        max_len_spacer = max(len(x) for x in self.list_spacers)
        padding_size_repeat = max_len_repeat_start + 3
        padding_size_spacer = max_len_spacer + 3

        for repeat_start, repeat, spacer, errors in zip(self.list_repeat_absolute_starts,
                                                        self.list_dotted_gaped_repeats,
                                                        self.list_spacers,
                                                        self.fuzzy_counts):
            repeat_start += 1
            line_padding_repeat = " " * (padding_size_repeat - len(str(repeat_start)))
            line_padding_spacer = " " * (padding_size_spacer - len(str(spacer)))
            substitutions, insertions, deletions = errors
            error_line = f"s:{substitutions} i:{insertions} d:{deletions}"
            line = f"{repeat_start}{line_padding_repeat}{repeat}   {spacer}{line_padding_spacer}{error_line}\n"
            dot_repr += line

        line_padding_last_repeat = " " * (padding_size_repeat - len(str(self.list_repeat_absolute_starts[-1])))
        last_repeat_start = self.list_repeat_absolute_starts[-1]
        last_repeat = self.list_dotted_gaped_repeats[-1]
        gap_instead_spacer = " " * padding_size_spacer
        l_substitutions, l_insertions, l_deletions = self.fuzzy_counts[-1]
        last_errors = f"s:{l_substitutions} i:{l_insertions} d:{l_deletions}"
        last_repeat_start += 1
        last_line = f"{last_repeat_start}{line_padding_last_repeat}{last_repeat}   {gap_instead_spacer}{last_errors}\n"
        dot_repr += last_line

        dot_repr += "_" * 100 + '\n'
        padding_consensus = " " * padding_size_repeat
        dot_repr += f"{padding_consensus}{self.gaped_matching_pattern}\n"

        dot_repr += "_" * 100 + '\n'
        padding_consensus_no_gap = " " * padding_size_repeat
        dot_repr += f"{padding_consensus_no_gap}{self.consensus}\n"

        dot_repr += "_" * 100 + '\n'
        padding_consensus_transitions = " " * padding_size_repeat
        dot_repr += f"{padding_consensus_transitions}{self.transition_consensus}\n"
        return dot_repr


class BlockRepeats:
    def __init__(self, consensus, repeats):
        self.consensus = consensus
        self.repeats = repeats

    def convert_block_into_columns(self):
        consesus_chars = list(self.consensus) + ["#"]
        columns_repeat = [list(column) for column in zip(*self.repeats)] + (["#"] * len(self.repeats))
        return [ColumnWithConsensusChar(con_char, l_repeats)
                for con_char, l_repeats in zip(consesus_chars, columns_repeat)]


class ColumnWithConsensusChar:
    def __init__(self, consensus_char, column_repeat_chars):
        self.consensus_char = consensus_char
        self.column_repeat_chars = column_repeat_chars

    def __repr__(self):
        column_repeat = "".join(self.column_repeat_chars)
        return f"{self.consensus_char}({column_repeat})"

    def __len__(self):
        return len(self.column_repeat_chars)


class BlockRecomputation:
    def __init__(self, blocks):
        self.blocks = blocks
        self.column_representations = [block.convert_block_into_columns() for block in self.blocks]
        self._recompute()

    def _insert_gaps(self, index):
        for c_r in self.column_representations:
            if c_r[index].consensus_char != " ":
                gap_column = ColumnWithConsensusChar(" ", [" " for _ in range(len(c_r[0]))])
                c_r.insert(index, gap_column)

    def _recompute(self):
        index = 0
        while index < len(self.column_representations[0]):
            if any([c_r[index].consensus_char == " " for c_r in self.column_representations]):
                self._insert_gaps(index)

            index += 1
        for c_r in self.column_representations:
            c_r.pop()

    def output_consensus_and_repeats(self):
        list_repeats = []
        for c_r in self.column_representations:
            for index in range(len(c_r[0])):
                repeat = [cwcc.column_repeat_chars[index] for cwcc in c_r]
                repeat = "".join(repeat)
                list_repeats.append(repeat)

        consensus = [cwcc.consensus_char for cwcc in self.column_representations[0]]
        consensus = "".join(consensus)

        return consensus, list_repeats

    def __repr__(self):
        return "\n".join(str(x) for x in self.column_representations)


def search_pair_handler(search_result, consensus_seq, relative_errors):
    list_consensus = list(consensus_seq)
    list_search_result = list(search_result)
    substitutions, insertions_orig, deletions = relative_errors

    insertions = apply_deletions_to_insertions(insertions_orig, deletions)

    insertions = set(insertions)
    #deletions = [d if d >= 0 else 0 for d in deletions] # normal solution
    if deletions:
        if min(deletions) < 0:
            deletions = [d - min(deletions) for d in deletions]
    deletions = set(deletions)
    same_index = insertions.intersection(deletions)

    insertions_deletions = sorted(list((insertions.union(deletions))))

    if same_index:
        insertions_deletions = set(insertions_deletions).difference(same_index)
        insertions_deletions = sorted(list(insertions_deletions))

    for _, ins_del in enumerate(insertions_deletions):
        if ins_del in insertions:
            list_consensus.insert(ins_del, " ")

        if ins_del in deletions:
            list_search_result.insert(ins_del, "-")

    consensus_gaped = "".join(list_consensus)
    search_result_with_deletions = "".join(list_search_result)

    return consensus_gaped, search_result_with_deletions


def apply_deletions_to_insertions(list_insertions, list_deletions):
    list_insertions_new = []
    for element in list_insertions:
        add = sum(element >= i for i in list_deletions)
        list_insertions_new.append(element + add)
    return list_insertions_new


def create_boundaries_for_intervals(intervals, flank_size):
    if intervals:
        start = max(0, intervals[0][0] - flank_size)
        end = intervals[0][1] + flank_size
        if len(intervals) == 1:
            return [(start, end)]
        else:
            ends = []
            starts = []
            for first, second in zip(intervals[:-1], intervals[1:]):
                end_first = first[1]
                begin_second = second[0]
                if end_first >= begin_second:
                    ends.append(end_first)
                    starts.append(begin_second)
                elif (end_first + 2 * flank_size) < begin_second:
                    ends.append(end_first + flank_size)
                    starts.append(begin_second - flank_size)
                else:
                    possible_flank = int((begin_second - end_first) / 2)
                    ends.append(end_first + possible_flank)
                    starts.append(begin_second - possible_flank)

            starts = [start] + starts
            ends.append(end)
            return [(st, end) for st, end in zip(starts, ends)]
    else:
        return []

#      Refine regex result (insertion instead of substitution)
#######################################################################
#######################################################################


class ArrayRefinerInsertionsDeletions:
    def __init__(self, crispr_candidate):
        self.crispr_candidate = crispr_candidate
        self.refined_candidate = crispr_candidate

        self._refine_insertions()

    def _refine_insertions(self):
        ari = ArrayRefinerInsertions(self.crispr_candidate)
        self.refined_candidate = ari.refined_candidate

    def _refine_deletions(self):
        pass

    def output(self):
        return self.refined_candidate


class ArrayRefinerInsertions:
    def __init__(self, crispr_candidate):
        self.crispr_candidate = crispr_candidate
        self.refined_candidate = crispr_candidate

        self.original_repeats = self.crispr_candidate.list_repeats
        self.original_repeats_gaped = self.crispr_candidate.list_repeats_gaped
        self.original_spacers = self.crispr_candidate.list_spacers
        self.original_starts = self.crispr_candidate.list_repeat_starts

        self.problematic_spacers_indexes = []

        self._search_potential_problems()
        self._refine_insertions()

    def _search_potential_problems(self):
        len_spacers = [len(s) for s in self.original_spacers]
        for index, len_spacer in enumerate(len_spacers):

            if index == 0:
                if len(len_spacers) > 2:
                    window_spacers_len = len_spacers[1:3]
                    if len(set(window_spacers_len)) == 1:
                        if window_spacers_len[0] != len_spacer:
                            self.problematic_spacers_indexes.append(index)
                else:
                    self.problematic_spacers_indexes.append(index)

            if index == 1:
                if len(len_spacers) > 2:
                    window_spacers_len = len_spacers[0:1] + len_spacers[2:3]
                    if len(set(window_spacers_len)) == 1:
                        if window_spacers_len[0] != len_spacer:
                            self.problematic_spacers_indexes.append(index)
                else:
                    self.problematic_spacers_indexes.append(index)

            if index > 1:
                if index not in self.problematic_spacers_indexes:
                    if index < (len(len_spacers) - 2):
                        window_spacers_len = len_spacers[index - 2:index] + len_spacers[index + 1:index + 3]
                        if len(set(window_spacers_len)) == 1:
                            if window_spacers_len[0] != len_spacer:
                                self.problematic_spacers_indexes.append(index)

            if index == (len(len_spacers) - 2):
                if index not in self.problematic_spacers_indexes:
                    if len(len_spacers) > 2:
                        window_spacers_len = len_spacers[len(len_spacers) - 3:len(len_spacers) - 2] + [len_spacers[-1]]
                        if len(set(window_spacers_len)) == 1:
                            if window_spacers_len[0] != len_spacer:
                                self.problematic_spacers_indexes.append(index)
                    else:
                        self.problematic_spacers_indexes.append(index)

            if index == (len(len_spacers) - 1):
                if index not in self.problematic_spacers_indexes:
                    if len(len_spacers) > 2:
                        window_spacers_len = len_spacers[len(len_spacers) - 2:len(len_spacers) - 1]
                        if len(set(window_spacers_len)) == 1:
                            if window_spacers_len[0] != len_spacer:
                                self.problematic_spacers_indexes.append(index)
                    else:
                        self.problematic_spacers_indexes.append(index)

    def _refine_insertions(self):
        for spacer_index in self.problematic_spacers_indexes:
            osir = OneSpacerInsertionRefiner(self.refined_candidate, spacer_index)
            self.refined_candidate = osir.refined_candidate


class RefineDeletions:
    def __init__(self, crispr_candidate):
        self.crispr_candidate = crispr_candidate

    def _refine_deletions(self):
        pass


class OneSpacerInsertionRefiner:
    def __init__(self, crispr_candidate, spacer_index):
        self.crispr_candidate = crispr_candidate
        self.spacer_index = spacer_index
        self.refined_candidate = crispr_candidate

        self.repeats = self.crispr_candidate.list_repeats
        self.repeats_gaped = self.crispr_candidate.list_repeats_gaped
        self.spacers = self.crispr_candidate.list_spacers
        self.starts = self.crispr_candidate.list_repeat_starts
        self.original_errors = self.crispr_candidate.total_mismatches

        self.indexes_insertions = None
        self.repeat_index_to_refine = None

        self._get_indexes_insertions()
        self._find_repeat_to_refine()

    def _get_indexes_insertions(self):
        self.indexes_insertions = list(set([index for repeat in self.repeats_gaped
                                            for index, char in enumerate(repeat) if char == ' ']))
        self.indexes_insertions = sorted(self.indexes_insertions)

    def _find_repeat_to_refine(self):
        self.left_repeat = self.repeats[self.spacer_index]
        self.right_repeat = self.repeats[self.spacer_index+1]

        if len(self.left_repeat) > len(self.right_repeat):
            self.repeat_index_to_refine = self.spacer_index
            self._refine_left_repeat()
        if len(self.left_repeat) < len(self.right_repeat):
            self.repeat_index_to_refine = self.spacer_index + 1
            self._refine_right_repeat()

    def _refine_left_repeat(self):
        gaped_repeat = self.repeats_gaped[self.repeat_index_to_refine]
        for global_insertion_index in self.indexes_insertions[::-1]:
            if gaped_repeat[global_insertion_index] != " ":
                index_last_insertion = global_insertion_index
                chars_in_the_gap = [r[index_last_insertion] for r in self.repeats_gaped]
                num_non_gaps = sum([1 if x != " " else 0 for x in chars_in_the_gap])
                insertion_only_in_the_refined_repeat = True if num_non_gaps == 1 else False

                index_last_insertion = global_insertion_index

                new_repeats = self.repeats[:]
                new_repeats[self.repeat_index_to_refine] = new_repeats[self.repeat_index_to_refine][:-1]

                new_spacers = self.spacers[:]
                new_spacers[self.spacer_index] = self.repeats[self.repeat_index_to_refine][-1] + \
                                                 new_spacers[self.repeat_index_to_refine]

                new_starts = self.starts[:]

                if insertion_only_in_the_refined_repeat:
                    new_repeats_gaped = []
                    for gaped_repeat in self.repeats_gaped:
                        if gaped_repeat[index_last_insertion] == " ":
                            gaped_repeat_new = gaped_repeat[:index_last_insertion] + \
                                               gaped_repeat[(index_last_insertion + 1):]
                            new_repeats_gaped.append(gaped_repeat_new)
                        else:
                            new_repeats_gaped.append(gaped_repeat)

                    new_repeats_gaped[self.repeat_index_to_refine] = new_repeats_gaped[self.repeat_index_to_refine][:-1]

                else:
                    new_repeats_gaped = self.repeats_gaped
                    repeat_to_refine = new_repeats_gaped[self.repeat_index_to_refine]
                    repeat_to_refine = repeat_to_refine[:index_last_insertion] + " " \
                                       + repeat_to_refine[index_last_insertion:-1]

                    new_repeats_gaped[self.repeat_index_to_refine] = repeat_to_refine

                new_candidate = CrisprCandidate(list_repeats=new_repeats,
                                                list_repeats_gaped=new_repeats_gaped,
                                                list_spacers=new_spacers,
                                                list_repeat_starts=new_starts)

                new_candidate_mismatches = new_candidate.total_mismatches
                if new_candidate_mismatches <= self.original_errors:
                    self.refined_candidate = new_candidate

                break

    def _refine_right_repeat(self):
        gaped_repeat = self.repeats_gaped[self.repeat_index_to_refine]
        for global_insertion_index in self.indexes_insertions:
            if gaped_repeat[global_insertion_index] != " ":
                index_first_insertion = global_insertion_index
                chars_in_the_gap = [r[index_first_insertion] for r in self.repeats_gaped]
                num_non_gaps = sum([1 if x != " " else 0 for x in chars_in_the_gap])
                insertion_only_in_the_refined_repeat = True if num_non_gaps == 1 else False

                new_repeats = self.repeats[:]
                new_repeats[self.repeat_index_to_refine] = new_repeats[self.repeat_index_to_refine][1:]

                new_spacers = self.spacers[:]
                new_spacers[self.spacer_index] = new_spacers[self.repeat_index_to_refine-1] + \
                                                 self.repeats[self.repeat_index_to_refine][0]

                new_starts = self.starts[:]
                new_starts[self.repeat_index_to_refine] = new_starts[self.repeat_index_to_refine] - 1

                if insertion_only_in_the_refined_repeat:
                    new_repeats_gaped = []
                    for gaped_repeat in self.repeats_gaped:
                        if gaped_repeat[index_first_insertion] == " ":
                            gaped_repeat_new = gaped_repeat[:index_first_insertion] + \
                                               gaped_repeat[(index_first_insertion + 1):]
                            new_repeats_gaped.append(gaped_repeat_new)
                        else:
                            new_repeats_gaped.append(gaped_repeat)

                    new_repeats_gaped[self.repeat_index_to_refine] = new_repeats_gaped[self.repeat_index_to_refine][1:]

                else:
                    new_repeats_gaped = self.repeats_gaped
                    repeat_to_refine = new_repeats_gaped[self.repeat_index_to_refine]
                    repeat_to_refine = repeat_to_refine[1:index_first_insertion+1] + " "\
                                       + repeat_to_refine[index_first_insertion+1:]
                    new_repeats_gaped[self.repeat_index_to_refine] = repeat_to_refine

                new_candidate = CrisprCandidate(list_repeats=new_repeats,
                                                list_repeats_gaped=new_repeats_gaped,
                                                list_spacers=new_spacers,
                                                list_repeat_starts=new_starts)

                new_candidate_mismatches = new_candidate.total_mismatches
                if new_candidate_mismatches <= self.original_errors:
                    self.refined_candidate = new_candidate

                break


class OneSpacerDeletionRefiner:
    def __init__(self, crispr_candidate, spacer_index):
        self.crispr_candidate = crispr_candidate
        self.spacer_index = spacer_index
        self.refined_candidate = crispr_candidate
