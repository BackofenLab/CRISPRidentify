import subprocess
import os
import re


class HMMMatch:
    def __init__(self, target, query, e_value, score):
        self.target = target
        self.query = query
        self.e_value = e_value
        self.score = score

    def __repr__(self):
        return "{} {} {} {}".format(self.target, self.query, self.e_value, self.score)


class HMMResultParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.list_matches = []
        self._input_parsing()

    def _input_parsing(self):
        with open(self.file_path, "r") as f:
            information_lines = f.readlines()

        information_lines = information_lines[3:]
        for line in information_lines:
            list_info = line.split()
            if list_info[0] != "#":
                score = float(list_info[5])
                e_value = float(list_info[4])
                target = list_info[0]
                query = list_info[2]

                hmm_match = HMMMatch(target, query, e_value, score)
                self.list_matches.append(hmm_match)

    def output(self):
        return self.list_matches


class CasGeneSearch:
    def __init__(self, file_name, hmm_model_file, e_value_th):
        self.file_name = file_name
        self.hmm_model_file = hmm_model_file
        self.e_value_th = e_value_th

        self.dict_protein_description = {}
        self.dict_protein_seq = {}
        self.hmm_matches = []
        self.dict_best_hmm_matches = {}
        self.dict_interval_cas = {}
        self.dict_cas_intervals = {}

        self._protein_extraction()
        self._hmm_extraction()
        self._e_value_filtration()
        self._final_result_creation()

    def _protein_extraction(self):
        print("Run protein extraction")
        cmd = 'tools/prodigal/prodigal'
        cmd += ' -i {} -p meta -a protein.fa'.format(self.file_name)

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        process.communicate()

        with open("protein.fa", "r") as f:
            lines = f.readlines()
        header_lines = [line[1:] for line in lines if ">" in line]
        keys = [line.split()[0] for line in header_lines]
        coordinates = [(int(line.split()[2]), int(line.split()[4])) for line in header_lines]
        strands = [int(line.split()[6]) for line in header_lines]
        strands = ["Forward" if x == 1 else "Reverse" for x in strands]
        values = [(*coordinate, strand) for coordinate, strand in zip(coordinates, strands)]

        self.dict_protein_description = {key: value for key, value in zip(keys, values)}

        header_lines_indexes = [index for index, line in enumerate(lines) if ">" in line]
        for index, header_index_first in enumerate(header_lines_indexes[:-1]):
            key = (int(lines[header_index_first].split()[2]), int(lines[header_index_first].split()[4]))
            header_index_second = header_lines_indexes[index+1]
            protein_seq = ""
            for internal_index in range(header_index_first+1, header_index_second):
                protein_seq += lines[internal_index].strip()

            self.dict_protein_seq[key] = protein_seq

        if header_lines_indexes:
            last_protein_st_end = (int(lines[header_lines_indexes[-1]].split()[2]),
                                   int(lines[header_lines_indexes[-1]].split()[4]))
            last_protein_seq = ""
            for internal_index in range(header_lines_indexes[-1], len(lines)):
                last_protein_seq += lines[internal_index].strip()
                self.dict_protein_seq[last_protein_st_end] = last_protein_seq

    def _hmm_extraction(self):
        print("Run hmm search")
        if os.stat('protein.fa').st_size != 0:
            cmd = 'tools/hmm_search/hmmsearch --tblout {} {} {}'.format('result_hmm.out', self.hmm_model_file,
                                                                        'protein.fa')
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            process.communicate()

            hmm_result_parser = HMMResultParser("result_hmm.out")
            self.hmm_matches = hmm_result_parser.output()
            os.remove('result_hmm.out')

        os.remove('protein.fa')

    def _e_value_filtration(self):
        print("Filtering by e value")
        self.hmm_matches = [hmm_match for hmm_match in self.hmm_matches if hmm_match.e_value <= self.e_value_th]

    def _final_result_creation(self):
        for hmm_match in self.hmm_matches:
            target = hmm_match.target
            query = hmm_match.query
            e_value = hmm_match.e_value
            if target not in self.dict_best_hmm_matches:
                self.dict_best_hmm_matches[target] = (e_value, query)
            else:
                current_e_val = self.dict_best_hmm_matches[target][0]
                if e_value < current_e_val:
                    self.dict_best_hmm_matches[target] = (e_value, query)

        for key, value in self.dict_best_hmm_matches.items():
            protein_description = self.dict_protein_description[key]
            protein_interval = protein_description[0], protein_description[1]
            strand = protein_description[2]
            e_value, query = value
            try:
                cas = re.search(r"cas\d+", query).group()
            except AttributeError:
                cas = query
                cas = cas.split("-")[0]

            self.dict_interval_cas[protein_interval] = cas, strand

        self.dict_cas_intervals = {cas: [key for key in self.dict_interval_cas.keys()
                                         if self.dict_interval_cas[key] == cas]
                                   for cas in self.dict_interval_cas.values()}

    def output_by_cas(self):
        return self.dict_cas_intervals

    def report_maker_by_cas(self, file_name):
        with open(file_name, "w") as f:
            for key, values in self.dict_cas_intervals.items():
                for value in values:
                    f.write("{}\t{}\t{}\t{}\n".format(key, value[0], value[1], self.dict_protein_seq[value]))

    def report_maker_by_interval(self, file_name):
        with open(file_name, "w") as f:
            for key in sorted(self.dict_interval_cas.keys()):
                value = self.dict_interval_cas[key]
                try:
                    cas = re.search(r"cas\d+", value).group()
                except AttributeError:
                    cas = value
                    cas = cas.split("-")[0]
                f.write("{}\t{}\n".format(key, cas))
