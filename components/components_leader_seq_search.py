from components_rev_com_crispr import rev_compliment_seq


class FullLeaderSeqSearch:
    def __init__(self, list_crispr_candidates, full_dna):
        self.list_crispr_candidates = list_crispr_candidates
        self.full_dna = full_dna

        self.dict_leader_data = {}

        self._compute_all_leaders()

    def _compute_all_leaders(self):
        self.dict_leader_data = {index: LeaderSeqSearch(crispr_candidate, self.full_dna).output()
                                 for index, crispr_candidate in enumerate(self.list_crispr_candidates)}

    def output(self):
        return self.dict_leader_data


class LeaderSeqSearch:
    def __init__(self, crispr_candidate, full_dna):
        self.crispr_candidate = crispr_candidate
        self.full_dna = full_dna

        self.leader = None

        self._compute_leader_seq()

    def _compute_leader_seq(self):
        list_repeat_indexes = self.crispr_candidate.list_repeat_starts
        list_repeats = self.crispr_candidate.list_repeats
        first_index = list_repeat_indexes[0]
        second_index = list_repeat_indexes[1]

        strand = "Forward" if first_index < second_index else "Reverse"

        if strand == "Forward":
            start = first_index - 101
            end = first_index - 1
            self.leader = self.full_dna[start:end+1]
        else:
            first_repeat = list_repeats[0]
            len_first_repeat = len(first_repeat)
            index_leader_start = first_index + len_first_repeat
            self.leader = self.full_dna[index_leader_start:index_leader_start+100]
            self.leader = rev_compliment_seq(self.leader)

    def output(self):
        return self.leader
