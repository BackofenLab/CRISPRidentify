import regex
from detection import FuzzySearch
from crispr_candidate import CrisprCandidate


class AdvancedFuzzySearch:
    def __init__(self, sequence, sequence_start, left_flank, left_flank_start,
                 right_flank, right_flank_start, repeat):

        self.pattern = "(?e)" + "({})".format(repeat) + "{i<=3,d<=3,s<=3,i+d+s<=6}"
        self.pattern_flank = "(?e)" + "({})".format(repeat) + "{i<=5,d<=5,s<=5,i+d+s<=10}"

        self.sequence = sequence
        self.sequence_start = sequence_start
        self.left_flank = left_flank
        self.left_flank_start = left_flank_start
        self.right_flank = right_flank
        self.right_flank_start = right_flank_start
        self.repeat = repeat
        self.repeat_candidate = repeat

        self._match()
        self._get_repeats()
        self._get_spacers()
        self._get_errors()
        self._get_absolute_starts()
        self._get_relative_error_indexes()
        self._get_all_insertions_deletions()
        self._get_gaped_repeats()
        self._compute_dot_representation()
        self.matching_object = None

    def _match(self):
        """Run fuzzy search on the sequence"""
        self.matching_object = list(regex.finditer(self.pattern, self.sequence, ENHANCEMATCH=True))
        self.matching_l_flank = list(regex.finditer(self.pattern_flank, self.left_flank, ENHANCEMATCH=True))
        self.matching_r_flank = list(regex.finditer(self.pattern_flank, self.right_flank, ENHANCEMATCH=True))

        self.matching_l_flank = self.matching_l_flank[0] if self.matching_l_flank else None
        self.matching_r_flank = self.matching_r_flank[-1] if self.matching_r_flank else None

        if self.matching_object:
            self.match_hit = True
        else:
            self.match_hit = False

    def _get_repeats(self):
        """Obtain the sequences that were matched"""
        self.list_repeats_main = [match.group() for match in self.matching_object]

        self.repeat_left = self.matching_l_flank.group() if self.matching_l_flank else None
        self.repeat_right = self.matching_r_flank.group() if self.matching_r_flank else None

        self.list_repeats = self.list_repeats_main
        if self.repeat_left:
            self.list_repeats = [self.repeat_left] + self.list_repeats
        if self.repeat_right:
            self.list_repeats = self.list_repeats + [self.repeat_right]

    def _get_spacers(self):
        """Obtain the spacers sequences"""
        spacer_starts = [match.end() for match in self.matching_object][:-1]
        spacer_end = [match.start() for match in self.matching_object][1:]
        spacer_intervals = [(start, end) for start, end in zip(spacer_starts, spacer_end)]
        self.list_spacers = [self.sequence[interval[0]:interval[1]] for interval in spacer_intervals]

        spacer_l_start = self.matching_l_flank.end() if self.matching_l_flank else None
        if spacer_l_start:
            spacer_l = self.left_flank[spacer_l_start:]
            self.list_spacers = [spacer_l] + self.list_spacers

        spacer_r_end = self.matching_r_flank.start() if self.matching_r_flank else None
        if spacer_r_end:
            spacer_r = self.right_flank[:spacer_r_end]
            self.list_spacers.append(spacer_r)

    def _get_errors(self):
        """Obtain all the errors"""
        self.number_substitutions = 0
        self.number_insertions = 0
        self.number_deletions = 0

        for match in self.matching_object:
            s, i, d = match.fuzzy_counts
            self.number_substitutions += s
            self.number_insertions += i
            self.number_deletions += d

        self.number_errors = self.number_substitutions + self.number_insertions + self.number_deletions

        self.errors_dot_repr = []
        for match in self.matching_object:
            self.errors_dot_repr.append(match.fuzzy_counts)

        if self.matching_l_flank:
            self.errors_dot_repr = [self.matching_l_flank.fuzzy_counts] + self.errors_dot_repr

        if self.matching_r_flank:
            self.errors_dot_repr = self.errors_dot_repr + [self.matching_r_flank.fuzzy_counts]

    def _get_absolute_starts(self):
        """Obtaining absolute start of matched sequences"""
        self.list_absolute_start = [self.sequence_start + match.start() for match in self.matching_object]

        if self.matching_l_flank:
            self.list_absolute_start = [self.left_flank_start + self.matching_l_flank.start()] +\
                                       self.list_absolute_start

        if self.matching_r_flank:
            self.list_absolute_start = self.list_absolute_start +\
                                       [self.right_flank_start + self.matching_r_flank.start()]

        self.start_end = '{}_{}'.format(self.list_absolute_start[0], self.list_absolute_start[-1])

    def _get_relative_error_indexes(self):
        """Obtaining indexes of errors which are related to each sequence and not global"""
        self.list_relative_error_indexes = []
        for match in self.matching_object:
            tuple_match_errors = match.fuzzy_changes
            list_relative_errors = [[e - match.start() for e in err_type] for err_type in tuple_match_errors]
            self.list_relative_error_indexes.append(list_relative_errors)

        if self.matching_l_flank:
            tuple_match_errors_l_flank = self.matching_l_flank.fuzzy_changes
            list_relative_errors = [[e - self.matching_l_flank.start() for e in err_type]
                                    for err_type in tuple_match_errors_l_flank]

            self.list_relative_error_indexes = [list_relative_errors] + self.list_relative_error_indexes

        if self.matching_r_flank:
            tuple_match_errors_r_flank = self.matching_r_flank.fuzzy_changes
            list_relative_errors = [[e - self.matching_r_flank.start() for e in err_type]
                                    for err_type in tuple_match_errors_r_flank]

            self.list_relative_error_indexes.append(list_relative_errors)

    def _get_all_insertions_deletions(self):
        self.list_all_insertions = []
        self.list_all_deletions = []

        for err_type in self.list_relative_error_indexes:
            self.list_all_insertions += self.unique_gaps_end(err_type[1])
            self.list_all_deletions += self.unique_gaps_end(err_type[2])

        self.list_all_insertions = list(set(self.list_all_insertions))
        self.list_all_deletions = list(set(self.list_all_deletions))

    @staticmethod
    def gap_introduction(sequence, list_of_gaps, gap_char):
        new_seq = ""
        for index in range(len(sequence) + len(list_of_gaps)):
            if index in list_of_gaps:
                new_seq += gap_char
            else:
                new_seq += sequence[0]
                sequence = sequence[1:]

        return new_seq

    @staticmethod
    def apply_insertions_to_deletions(list_insertions, list_deletions):
        list_deletions_new = []
        for element in list_deletions:
            add = sum(element >= i for i in list_insertions)
            list_deletions_new.append(element + add)
        return list_deletions_new

    @staticmethod
    def unique_gaps_end(list_gaps):
        list_gaps_new = []
        for el in list_gaps:
            while True:
                if el in list_gaps_new:
                    el += 1
                else:
                    list_gaps_new.append(el)
                    break
        return list_gaps_new

    @staticmethod
    def regex_bug_fix(list_gaps):
        return [x if x >= 0 else 0 for x in list_gaps]

    def _get_gaped_repeats(self):
        self.list_gaped_repeats = []
        for index, repeat in enumerate(self.list_repeats):
            local_insertions = self.list_relative_error_indexes[index][1]
            local_insertions = self.regex_bug_fix(local_insertions)
            local_insertions = self.unique_gaps_end(local_insertions)
            insertions = list(set(self.list_all_insertions).difference(set(local_insertions)))

            deletions = self.list_relative_error_indexes[index][2]
            deletions = self.regex_bug_fix(deletions)
            deletions = self.unique_gaps_end(deletions)

            gaped_repeat = self.gap_introduction(repeat, deletions, "-")
            gaped_repeat = self.gap_introduction(gaped_repeat, insertions, " ")

            self.list_gaped_repeats.append(gaped_repeat)
        self.gaped_matching_pattern = self.gap_introduction(self.repeat_candidate, self.list_all_insertions, " ")

    def __repr__(self):
        return self.dot_representation

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

    def _compute_dot_representation(self):
        string = ""
        if self.matching_object:
            max_length_start_index = max(len(str(start)) for start in self.list_absolute_start) + 3
            try:
                max_length_spacer = max(len(spacer) for spacer in self.list_spacers) + 3
            except ValueError:
                max_length_spacer = 0

            for index, gaped_repeat in enumerate(self.list_gaped_repeats):
                repeat_start_index = self.list_absolute_start[index] + 1
                n_gaps_after_start = max_length_start_index - len(str(self.list_absolute_start[index]))

                if index < len(self.list_spacers):
                    spacer = self.list_spacers[index]
                else:
                    spacer = ""
                n_gaps_after_spacer = max_length_spacer - len(spacer)

                s, i, d = self.errors_dot_repr[index]
                errors = "s:{} i:{} d:{}".format(s, i, d)
                dotted_repeats = self.apply_dots(gaped_repeat, self.gaped_matching_pattern)

                string += "{}{}{}  {}{}{}\n".format(repeat_start_index,
                                                      " " * n_gaps_after_start,
                                                      dotted_repeats,
                                                      spacer,
                                                      " " * n_gaps_after_spacer,
                                                      errors)

            string += "_" * 100 + "\n"

            string += " " * max_length_start_index + self.gaped_matching_pattern
            string += " " * (max_length_spacer + 2) + "s:{} i:{} d:{}".format(self.number_substitutions,
                                                                              self.number_insertions,
                                                                              self.number_deletions) + "\n"

            #string += "_" * 100 + "\n"
            #string += " " * max_length_start_index + self.repeat_candidate

        self.dot_representation = string

    def dot_repr(self):
        return self.dot_representation


class MockArray:
    def __init__(self, start, consensus):
        self.list_repeat_starts = [start]
        self.consensus = consensus


class FullDegeneratedSearch:
    def __init__(self, full_dna, list_crisprs):
        self.full_dna = full_dna
        self.list_crisprs = list_crisprs

        self.dict_degenerated = {}
        self._compute_all_degenerated()

    def _compute_all_degenerated(self):
        for index, crispr_array in enumerate(self.list_crisprs):
            ds = DegeneratedSearch(self.full_dna, crispr_array)
            self.dict_degenerated[index] = ds.output()

    def output(self):
        return self.dict_degenerated


class DegeneratedSearch:
    def __init__(self, full_dna, crispr_array):
        self.full_dna = full_dna[0:int(len(full_dna)/2)]
        self.crispr_array = crispr_array

        self.degenerate_results = {}

        self._get_flank()
        self._compute_degenerated()

    def _get_flank(self):
        self.crispr_start = self.crispr_array.list_repeat_starts[0]
        self.crispr_end = self.crispr_array.list_repeat_starts[-1] + len(self.crispr_array.list_repeats[-1])
        self.crispr_consensus = self.crispr_array.consensus

        self.l_flank_start = max(0, self.crispr_start - 150)
        self.r_flank_end = min(len(self.full_dna) - 1, self.crispr_end + 150)
        self.l_flank = self.full_dna[self.l_flank_start:self.crispr_start]
        self.r_flank = self.full_dna[self.crispr_end:self.r_flank_end+1]
        self.sequence = self.full_dna[self.crispr_start:self.crispr_end]

    def _compute_degenerated(self):
        fuzzy = AdvancedFuzzySearch(sequence=self.sequence, sequence_start=self.crispr_start,
                                    left_flank=self.l_flank, left_flank_start=self.l_flank_start,
                                    right_flank=self.r_flank, right_flank_start=self.crispr_end,
                                    repeat=self.crispr_consensus)

        self.crispr_candidate = CrisprCandidate(fuzzy.list_repeats, fuzzy.list_gaped_repeats,
                                                fuzzy.list_spacers, fuzzy.list_absolute_start)

        self.dot_representation = fuzzy.dot_repr()

    def output(self):
        return self.crispr_candidate

    def output_dot_repr(self):
        return self.dot_representation


