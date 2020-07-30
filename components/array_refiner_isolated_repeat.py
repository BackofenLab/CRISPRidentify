from crispr_candidate import CrisprCandidate


class ArrayRefinerIsolatedRepeat:
    def __init__(self, crispr_candidate, spacer_len_th):
        self.crispr_candidate = crispr_candidate
        self.spacer_len_th = spacer_len_th

        self._remove_isolated_recursive()

    def _remove_isolated_recursive(self):

        while 1:
            list_spacers = self.crispr_candidate.list_spacers
            index_long_spacers = [index for index, spacer in enumerate(list_spacers)
                                  if len(spacer) > self.spacer_len_th]

            list_repeats = self.crispr_candidate.list_repeats
            list_repeats_gaped = self.crispr_candidate.list_repeats_gaped
            list_repeat_starts = self.crispr_candidate.list_repeat_starts

            if len(list_repeats) < 3:
                return

            if 0 in index_long_spacers:
                list_repeats = list_repeats[1:]
                list_spacers = list_spacers[1:]
                list_repeats_gaped = list_repeats_gaped[1:]
                list_repeat_starts = list_repeat_starts[1:]
                self.crispr_candidate = CrisprCandidate(list_repeats, list_repeats_gaped,
                                                        list_spacers, list_repeat_starts)

            elif (len(list_spacers) - 1) in index_long_spacers:
                list_repeats = list_repeats[:-1]
                list_spacers = list_spacers[:-1]
                list_repeats_gaped = list_repeats_gaped[:-1]
                list_repeat_starts = list_repeat_starts[:-1]
                self.crispr_candidate = CrisprCandidate(list_repeats, list_repeats_gaped,
                                                        list_spacers, list_repeat_starts)
            else:
                return

    def output(self):
        return self.crispr_candidate
