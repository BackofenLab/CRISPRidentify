from components_strand import get_orientation


class ResultParser:
    def __init__(self, file_name):
        self.file_name = file_name
        self._parse_results()

    def _parse_results(self):
        with open(self.file_name, "r") as f:
            lines = f.readlines()
        headers = [line for line in lines if "number of Repeats" in line]

        categories = [header.split(":")[0].split(" ")[0] for header in headers]

        crispr_indexes = [header.split(",")[0].split(":")[1][1:] for header in headers]

        start_ends = [header.split(",")[1][1:] for header in headers]
        starts = [start_end.split("-")[0] for start_end in start_ends]
        ends = [start_end.split("-")[1] for start_end in start_ends]

        numbers_of_repeats = [header.split(",")[2].split(" ")[-1] for header in headers]

        avg_spacer_lengths = [header.split(":")[-1][1:].strip() for header in headers]

        indexes_consensus = [index + 1 for index, line in enumerate(lines) if "_______" in line]
        consensus_repeats = [lines[index].split("s")[0].strip().replace(" ", "") for index in indexes_consensus]

        self.list_results = [CRISPRArrayParsingContainer(category, crispr_index, start, end,
                                                         number_of_repeats, avg_spacer_length, consensus_repeats)
                             for category, crispr_index, start, end,
                             number_of_repeats, avg_spacer_length, consensus_repeats in
                             zip(categories, crispr_indexes, starts, ends,
                             numbers_of_repeats, avg_spacer_lengths, consensus_repeats)]

    def output(self):
        return self.list_results


class CRISPRArrayParsingContainer:
    def __init__(self, category, crispr_index, start, end, number_of_repeats, avg_spacer_length, consensus_repeat):
        self.category = category if category != "CRISPR" else "Best"
        self.crispr_index = crispr_index
        self.start = start
        self.end = end
        self.number_of_repeats = number_of_repeats
        self.avg_spacer_length = avg_spacer_length
        self.consensus_repeat = consensus_repeat

        self.strand = get_orientation(self.consensus_repeat)
        self.strand = "Forward" if self.strand else "Reverse"

        if self.strand == "Reverse":
            self.start, self.end = self.end, self.start

    def __repr__(self):
        return f"{self.category}-{self.start}-{self.end}"
