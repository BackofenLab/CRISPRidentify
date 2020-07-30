import regex
import sys


from crispr_candidate import CrisprCandidate


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
            if index > 1:
                if index < (len(len_spacers) - 3):
                    window_spacers_len = len_spacers[index - 2:index] + len_spacers[index + 1:index + 3]
                    if len(set(window_spacers_len)) == 1:
                        if window_spacers_len[0] != len_spacer:
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

                new_repeats = self.repeats[:]
                new_repeats[self.repeat_index_to_refine] = new_repeats[self.repeat_index_to_refine][:-1]

                new_spacers = self.spacers[:]
                new_spacers[self.spacer_index] = new_spacers[self.repeat_index_to_refine] + \
                                                 new_repeats[self.repeat_index_to_refine][-1]

                new_starts = self.starts[:]

                new_repeats_gaped = []
                for gaped_repeat in self.repeats_gaped:
                    if gaped_repeat[index_last_insertion] == " ":
                        gaped_repeat_new = gaped_repeat[:index_last_insertion] + \
                                           gaped_repeat[(index_last_insertion + 1):]
                        new_repeats_gaped.append(gaped_repeat_new)
                    else:
                        new_repeats_gaped.append(gaped_repeat)

                new_repeats_gaped[self.repeat_index_to_refine] = new_repeats_gaped[self.repeat_index_to_refine][:-1]

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

                new_repeats = self.repeats[:]
                new_repeats[self.repeat_index_to_refine] = new_repeats[self.repeat_index_to_refine][1:]

                new_spacers = self.spacers[:]
                new_spacers[self.spacer_index] = new_repeats[self.repeat_index_to_refine][0] + \
                                                 new_spacers[self.repeat_index_to_refine]

                new_starts = self.starts[:]
                new_starts[self.repeat_index_to_refine] = new_starts[self.repeat_index_to_refine] - 1

                new_repeats_gaped = []
                for gaped_repeat in self.repeats_gaped:
                    if gaped_repeat[index_first_insertion] == " ":
                        gaped_repeat_new = gaped_repeat[:index_first_insertion] + \
                                           gaped_repeat[(index_first_insertion + 1):]
                        new_repeats_gaped.append(gaped_repeat_new)
                    else:
                        new_repeats_gaped.append(gaped_repeat)

                new_repeats_gaped[self.repeat_index_to_refine] = new_repeats_gaped[self.repeat_index_to_refine][1:]

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



