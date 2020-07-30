import collections
import json


class CrisprConsensus(object):
    def __init__(self, list_repeats_gaped):
        self.list_repeats_gaped = list_repeats_gaped        
        
        self.num_different_repeat_length = None
        self.consensus = None
        self.consensus_no_gap = None
        self.len_consensus = None
        self.number_repeats = None
        
        self._check_repeat_length()        
        self._compute_consensus()        
    
    def _check_repeat_length(self):
        list_lengths = [len(repeat) for repeat in self.list_repeats_gaped]
        self.num_different_repeat_length = len(set(list_lengths))        
    
    def _compute_consensus(self):
        if self.num_different_repeat_length == 0:        
            print('Got repeats of 0 length')
        elif self.num_different_repeat_length != 1:
            print('Got a case with different repeat lengths')
        else:            
            self.consensus = ''
            for char_ind, _ in enumerate(self.list_repeats_gaped[0]):
                list_char_in_column = [repeat[char_ind] for repeat in self.list_repeats_gaped]
                counter = collections.Counter(list_char_in_column)
                freq = counter.most_common()                
                most_common_char = freq[0][0] if freq[0][0] != '-' else freq[1][0]
                self.consensus += most_common_char

        self.consensus_no_gap = self.consensus.replace(' ', '').replace('+', '')
        self.len_consensus = len(self.consensus_no_gap)
    
    def output(self):
        return self.consensus_no_gap, self.consensus


class CrisprCandidate(object):
    def __init__(self, list_repeats, list_repeats_gaped, list_spacers, list_repeat_starts):
        self.list_repeats = list_repeats
        self.list_repeats_gaped = list_repeats_gaped
        self.list_spacers = list_spacers
        self.list_repeat_starts = list_repeat_starts
        
        self.list_repeat_mismatches = []
        self.list_mismatches_indexes = []
                
        self.consensus = None
        self.consensus_gaped = None
        self.total_mismatches = None
        
        self._filter_redundant_insertion_deletions()        
        self._compute_consensus()
        self._compute_mismatches()
        
    def _filter_redundant_insertion_deletions(self):
        def _fix_repeats(list_repeats, list_bad_indexes_to_fix):
            list_repeats_new = []
            for repeat in list_repeats:
                list_repeats_new.append(_fix_repeat(repeat, list_bad_indexes_to_fix))
                
            return list_repeats_new        
        
        def _fix_repeat(repeat, list_bad_indexes_to_fix):
            new_repeat = ''
            for index, char in enumerate(repeat):
                if index not in list_bad_indexes_to_fix:
                    new_repeat += char
            
            return new_repeat
            
        list_bad_indexes = []
        for char_ind, _ in enumerate(self.list_repeats_gaped[0]):            
            list_char_in_column = [repeat[char_ind] for repeat in self.list_repeats_gaped]
            chars = set(list_char_in_column)
            
            if chars == {' '} or chars == {'-'}:
                list_bad_indexes.append(char_ind)       
                
        if list_bad_indexes:           
            self.list_repeats_gaped = _fix_repeats(self.list_repeats_gaped, list_bad_indexes)            
        
    def _compute_consensus(self):
        self.consensus, self.consensus_gaped = CrisprConsensus(self.list_repeats_gaped).output()
        
    def _compute_mismatches(self):
        def _compute_mismatches_repeat(gaped_repeat):
            substitutions = 0
            insertions = 0
            deletions = 0
            list_mismatches_indexes_one_repeat = []
            for index, char_repeat, char_con_repeat in zip(range(len(gaped_repeat)),
                                                           gaped_repeat,
                                                           self.consensus_gaped):
                
                if char_con_repeat == ' ':                
                    if char_repeat != ' ':
                        insertions += 1
                        list_mismatches_indexes_one_repeat.append(index)
                else:
                    if char_repeat == char_con_repeat:
                        pass
                    else:                        
                        if char_repeat == '-':
                            deletions += 1
                            list_mismatches_indexes_one_repeat.append(index)
                        elif char_repeat == ' ':
                            deletions += 1
                        else:
                            substitutions += 1
                            list_mismatches_indexes_one_repeat.append(index)
                             
            return substitutions, insertions, deletions, list_mismatches_indexes_one_repeat
            
        for gaped_repeat in self.list_repeats_gaped:
            s, i, d, list_mismatches_indexes_one_repeat = _compute_mismatches_repeat(gaped_repeat)
            total = s + i + d
            repeat_stats = [s, i, d, total]
            self.list_repeat_mismatches.append(repeat_stats)
            self.list_mismatches_indexes.append(list_mismatches_indexes_one_repeat)

        self.total_mismatches = sum([x[3] for x in self.list_repeat_mismatches])
    
    def dot_repeat(self, gaped_repeat):
        string = ''
        substitutions = 0
        insertions = 0
        deletions = 0
        for char_repeat, char_consensus in zip(gaped_repeat, self.consensus_gaped):
            if char_consensus == ' ':                
                string += char_repeat
                if char_repeat != ' ':
                    insertions += 1
            else:
                if char_repeat == char_consensus:
                    string += '.'
                else:
                    string += char_repeat
                    if char_repeat == '-':
                        deletions += 1
                    elif char_repeat == ' ':
                        deletions += 1
                    else:
                        substitutions += 1
        return string, substitutions, insertions, deletions
        
    def dot_repr(self):
        string = ''
        g_s, g_i, g_d = 0, 0, 0
        max_length_start_index = max(len(str(start)) for start in self.list_repeat_starts) + 3
        max_length_spacer = max(len(spacer) for spacer in self.list_spacers) + 3
        
        for index, gaped_repeat in enumerate(self.list_repeats_gaped):
            repeat_start_index = self.list_repeat_starts[index] + 1
            n_gaps_after_start = max_length_start_index - len(str(repeat_start_index))           

            if index == len(self.list_spacers):
                spacer = ""
            else:
                spacer = self.list_spacers[index]
            n_gaps_after_spacer = max_length_spacer - len(spacer)            
            
            dotted_repeats, s, i, d = self.dot_repeat(gaped_repeat)
            errors = "   s:{} i:{} d:{}".format(s, i, d)
            g_s += s
            g_i += i
            g_d += d

            string += "{}{}{}  {}{}{}\n".format(repeat_start_index,
                                                " " * n_gaps_after_start,
                                                dotted_repeats, spacer,
                                                " " * n_gaps_after_spacer,
                                                errors)

        string += "_" * 100 + "\n"       

        string += " " * max_length_start_index + self.consensus_gaped
        string += " " * (max_length_spacer + 2) + "   s:{} i:{} d:{}".format(g_s, g_i, g_d) + "\n"
        
        return string

    def write_file(self, file_name):
        with open(file_name, "w") as f:
            f.write(self.dot_repr())

    def write_as_json(self, filename):
        dict_to_write = {"repeat_begins": self.list_repeat_starts,
                         "repeats": self.list_repeats,
                         "repeats_gaped": self.list_repeats_gaped,
                         "spacers": self.list_spacers}

        with open(filename, 'w') as outfile:
            json.dump(dict_to_write, outfile)

    def compute_stats(self):
        start = self.list_repeat_starts[0] + 1
        end = self.list_repeat_starts[-1] + len(self.list_repeats[-1])
        avg_repeat = len(self.consensus)
        avg_spacer = int(sum((len(spacer) for spacer in self.list_spacers)) / len(self.list_spacers))
        number_repeats = len(self.list_repeats)
        return {"start": start, "end": end, "avg_repeat": avg_repeat,
                "avg_spacer": avg_spacer, "number_repeats": number_repeats}

    @classmethod
    def init_from_json(cls, file_name):
        with open(file_name) as json_file:
            dict_data = json.load(json_file)

            list_repeas = dict_data["repeats"]
            list_repeats_starts = dict_data["repeat_begins"]
            list_spacers = dict_data["spacers"]
            list_repeats_gaped = dict_data["repeats_gaped"]

        return cls(list_repeats=list_repeas, list_spacers=list_spacers,
                   list_repeats_gaped=list_repeats_gaped, list_repeat_starts=list_repeats_starts)
