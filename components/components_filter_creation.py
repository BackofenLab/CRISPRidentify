import os
import subprocess


class FilterApproximation:
    def __init__(self, list_repeats):
        self.list_repeats = list_repeats
        self.list_repeats = sorted(self.list_repeats)
        self.max_seq = None
        self.min_seq = None
        self.list_missing_candidates = []

        self._relative_path_generation()
        self._make_fasta_with_repeats()

        if len(self.list_repeats) > 2:
            self._obtain_max_sequence()
            self._obtain_min_sequence()
            self._compute_all_the_missing_cases()
            self._report_results()

        self._clean_up()

    def _relative_path_generation(self):
        full_path = os.path.realpath(__file__)
        self.absolute_path_to_tools = "/" + "/".join(full_path.split("/")[:-2]) + "/tools/"

    def _make_fasta_with_repeats(self):
        with open("clustal_repeats.fa", "w") as f:
            lines = "".join([f">r{index}\n{repeat}\n" for index, repeat in enumerate(self.list_repeats, 1)])
            f.write(lines)

    def _run_clustal_omega_repeats(self):
        cmd = self.absolute_path_to_tools + "clustalOmega/clustalo -i clustal_repeats.fa -o loc_align.txt"
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

    def _find_max_clustal_omega_repeat_results(self):
        with open("loc_align.txt") as f:
            self.align_result_repeats = f.readlines()
        if self.align_result_repeats:
            if set([line.strip()[0] for line in self.align_result_repeats[::2]]) == {">"}:
                alignments = [line.strip() for line in self.align_result_repeats[1::2]]
                columns = [[al[index] for al in alignments] for index in range(len(alignments[0]))]
                columns = [[x for x in column if x != "-"]for column in columns]
                columns = [sorted(column) for column in columns]
                #most_common = [max(set(lst).difference(set("-")), key=lst.count) for lst in columns]
                most_common = [max(lst, key=lst.count) for lst in columns]
                self.max_seq = "".join(most_common)
        os.remove("loc_align.txt")

    def _find_min_clustal_omega_repeat_results(self):
        index_start = None
        index_end = None

        with open("loc_align.txt") as f:
            self.align_result_repeats = f.readlines()
        if self.align_result_repeats:
            if set([line.strip()[0] for line in self.align_result_repeats[::2]]) == {">"}:
                alignments = [line.strip() for line in self.align_result_repeats[1::2]]
                columns = [[al[index] for al in alignments] for index in range(len(alignments[0]))]
                columns = [[x for x in column if x != "-"] for column in columns]
                columns = [sorted(column) for column in columns]
                most_common = [max(lst, key=lst.count) for lst in columns]
                conserved = "".join(most_common)

                for index in range(len(alignments[0])):
                    column = [al[index] for al in alignments]
                    if "-" not in column:
                        if not index_start:
                            index_start = index
                        index_end = index
                if index_start and index_end:
                    self.min_seq = conserved[index_start:index_end+1]
        os.remove("loc_align.txt")

    def _run_clustal_omega_min_max(self):
        with open("min_max.fa", "w") as f:
            f.write(">max\n")
            f.write(f"{self.max_seq}\n")
            f.write(">min\n")
            f.write(f"{self.min_seq}\n")

        cmd = self.absolute_path_to_tools + "clustalOmega/clustalo -i min_max.fa -o min_max_align.txt"
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

    def _obtain_max_sequence(self):
        self._run_clustal_omega_repeats()
        self._find_max_clustal_omega_repeat_results()

    def _obtain_min_sequence(self):
        self._run_clustal_omega_repeats()
        self._find_min_clustal_omega_repeat_results()

    def _compute_all_the_missing_cases(self):
        if self.max_seq and self.min_seq:
            if len(self.min_seq) > 4:
                if abs(len(self.max_seq) - len(self.min_seq)) < 16:
                    length_longest = len(self.max_seq)
                    length_shortest = len(self.min_seq)
                    all_possible_substrings = [self.max_seq[i:j + 1] for i in range(length_longest)
                                               for j in range(i, length_longest)
                                               if len(self.max_seq[i:j + 1]) > length_shortest]

                    for candidate in all_possible_substrings:
                        if self.min_seq in candidate:
                            self.list_missing_candidates.append(candidate)

    def _report_results(self):
        missing = ""
        if len(self.list_missing_candidates) > 0:
            with open("missing.fa", "w") as f:
                for index, seq in enumerate(self.list_missing_candidates, 1):
                    f.write(f">{index}\n")
                    f.write(f"{seq}\n")

            if len(self.list_missing_candidates) > 1:
                cmd = self.absolute_path_to_tools + "clustalOmega/clustalo -i missing.fa -o missing_align.txt"
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                process.communicate()

                with open("missing_align.txt") as fr:
                    missing = fr.readlines()
                    missing = missing[1::2]

                if not missing:
                    missing = self.list_missing_candidates[0]

                os.remove("missing_align.txt")
            os.remove("missing.fa")

        """
        with open("log_filter.txt", "a") as fw:
            fw.write("Vmatch_results:\n")
            for repeat in self.list_repeats:
                fw.write(f"{repeat}\n")

            fw.write("\n")
            fw.write("RepeatAlignment:\n")
            fw.writelines(self.align_result_repeats)
            fw.write("\n")
            fw.write("Max_seq:\n")
            fw.write(f"{str(self.max_seq)}\n")
            fw.write("Min_seq:\n")
            fw.write(f"{str(self.min_seq)}\n")
            fw.write("Missing repeats:\n")
            fw.writelines(missing)
            fw.write("\n\n")
            fw.write("="*30)
            fw.write("\n")
            
        """

    def _clean_up(self):
        try:
            os.remove("clustal_repeats.fa")
        except FileNotFoundError:
            pass

    def output(self):
        return self.list_missing_candidates


if __name__ == "__main__":
    list_repeats_test = ["GAATCTTCGAGATAGAATTGCAAG", "CTTCGAGATAGAATTGCAAGGAAT", "CCCTTGAATCTTCGAGATAGAATTGCAAG"]
    fa = FilterApproximation(list_repeats_test)