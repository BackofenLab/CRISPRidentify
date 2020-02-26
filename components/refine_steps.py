from crispr_candidate import CrisprCandidate


class PartSpacerRepeatChecker:
    def __init__(self, crispr_candidate):
        self.crispr_candidate = crispr_candidate

    def _need_modification_check(self, min_interval_len):
        seq_to_check = self.crispr_candidate.consensus[:min_interval_len + 1]
        spacers = self.crispr_candidate.list_spacers
        for spacer in spacers:
            if seq_to_check in spacer:
                return True
        return False

    def _refine_array(self, min_interval_len):

        def split_spacer(spacer, seed):
            consensus = self.crispr_candidate.consensus
            if seed not in spacer:
                return spacer, "", ""
            else:
                seed_len = len(seed)
                while consensus[:seed_len+1] in spacer:
                    seed_len += 1

                index_of_seed = spacer.find(consensus[:seed_len])
                first_spacer = spacer[:index_of_seed]
                short_repeat = consensus[:seed_len]
                second_spacer = spacer[index_of_seed+seed_len:]

            return first_spacer, short_repeat, second_spacer

        def gap_short_repeat(short_repeat):
            gaps_indexes = [index for index, char in enumerate(old_consensus_gaped) if char == " "]
            gaped_short_repeat = ""
            index_on_short = 0

            for index, _ in enumerate(old_consensus):
                if index in gaps_indexes:
                    gaped_short_repeat += " "
                else:
                    gaped_short_repeat += short_repeat[index_on_short]
                    index_on_short += 1

                if index_on_short == len(short_repeat):
                    break
            gaped_short_repeat += " " * (len(old_consensus_gaped) - len(gaped_short_repeat))

            return gaped_short_repeat

        old_list_repeats_gaped = self.crispr_candidate.list_repeats_gaped
        old_list_repeats = self.crispr_candidate.list_repeats
        old_list_spacers = self.crispr_candidate.list_spacers
        old_list_repeats_starts = self.crispr_candidate.list_repeat_starts
        old_consensus_gaped = self.crispr_candidate.consensus_gaped
        old_consensus = self.crispr_candidate.consensus

        seed = old_consensus[:min_interval_len]
        list_new_pairs = [split_spacer(spacer, seed) for spacer in old_list_spacers]

        new_list_repeats_gaped = []
        new_list_repeats = []
        new_list_spacers = []
        new_list_repeats_starts = []

        for index, trio in enumerate(list_new_pairs):
            new_list_repeats_gaped.append(old_list_repeats_gaped[index])
            new_list_repeats.append(old_list_repeats[index])
            new_list_repeats_starts.append(old_list_repeats_starts[index])

            if not trio[1]:
                new_list_spacers.append(trio[0])
            else:
                new_list_spacers.append(trio[0])

                new_list_repeats.append(trio[1])
                new_list_repeats_gaped.append(gap_short_repeat(trio[1]))
                new_list_repeats_starts.append(old_list_repeats_starts[index] + len(trio[0]) + len(new_list_repeats[-2]))

                new_list_spacers.append(trio[2])

        return CrisprCandidate(list_repeats=new_list_repeats, list_repeats_gaped=new_list_repeats_gaped,
                               list_repeat_starts=new_list_repeats_starts, list_spacers=new_list_spacers)

    def naive_modification(self, min_interval_len):
        if self._need_modification_check(min_interval_len):
            return self._refine_array(min_interval_len)
        else:
            return None


class RepeatExtension:
    def __init__(self, crispr_candidate, percentage_identity, full_dna):
        self.crispr_candidate = crispr_candidate
        self.percentage_identity = percentage_identity
        self.full_dna = full_dna

    def extend_to_the_left(self):
        pass

    def extend_to_the_right(self):
        pass


class Refine:
    def search_for_part_spacers_repeat(self, list_crispr_candidates, min_repeat):
        list_found_cases = []
        for candidate in list_crispr_candidates:
            new_crispr = PartSpacerRepeatChecker(candidate).naive_modification(min_interval_len=min_repeat)
            if new_crispr:
                list_found_cases.append(new_crispr)
        return list_found_cases

    def search_for_repeat_extensions(self, list_crispr_candidates, number):
        pass


