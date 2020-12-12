import collections
import json
from functools import wraps
from itertools import groupby


class SameStartEndFilter:
    def __init__(self, dict_crispr_candidates):
        self.dict_crispr_candidates = dict_crispr_candidates
        self.dict_filtered_start_end_crispr_candidates = {}

        self._filter_fuzzy_searches_same_start_end()

    def _filter_fuzzy_searches_same_start_end(self):
        for cluster_seq, list_fuzzy_s in self.dict_crispr_candidates.items():
            list_start_end = [fuzzy_s.start_end for fuzzy_s in list_fuzzy_s]
            pattern_len = [len(fuzzy_s.repeat_candidate) for fuzzy_s in list_fuzzy_s]
            tuples_st_end_len = zip(list_start_end, pattern_len)

            list_categories = [[fuzzy_s for fuzzy_s in list_fuzzy_s if
                                (fuzzy_s.start_end, len(fuzzy_s.repeat_candidate)) == tuple_info]
                               for tuple_info in tuples_st_end_len]

            best_fuzzy_s = [sorted(category, key=lambda x: x.number_errors)[0]
                            for category in list_categories]

            best_fuzzy_s_unique_repeat = []
            u_repeats = []
            for b_fuz in best_fuzzy_s:
                repeat = b_fuz.repeat_candidate
                if repeat not in u_repeats:
                    u_repeats.append(repeat)
                    best_fuzzy_s_unique_repeat.append(b_fuz)

            self.dict_filtered_start_end_crispr_candidates[cluster_seq] = best_fuzzy_s_unique_repeat

    def output(self):
        return self.dict_filtered_start_end_crispr_candidates


#          For filtering out non CRISPR cases
#############################################################
#############################################################
DEBUG_MODE = False


def exception_handler(function):
    @wraps(function)
    def wrapper(*args, **kwargs):
        try:
            result = function(*args, **kwargs)
            return result
        except Exception:
            return False
    return wrapper


def printing_if_filtered(function):
    @wraps(function)
    def wrapper(*args, **kwargs):
        result = function(*args, **kwargs)
        if DEBUG_MODE:
            if not result:
                with open("filtered_results.txt", "a") as f:
                    f.write("\n\n")
                    f.write("\n".join([str(arg) for arg in args]))
                    f.write("\n\n")
                    f.write(function.__name__)
                    f.write("\n\n")

        return result
    return wrapper


class AdvancedFuzzySearchFilter:
    def __init__(self, min_column_dominance_repeat, min_avg_spacer_length,
                 max_spacer_length, max_column_dominance_spacer, max_allowed_consecutive_spacers,
                 max_allowed_same_spacers, max_inconsistent_columns, min_avg_repeat_length,
                 max_avg_repeat_length, max_avg_spacer_length, min_repeats):

        self.column_dominance = min_column_dominance_repeat
        self.min_avg_spacer_length = min_avg_spacer_length
        self.max_spacer_length = max_spacer_length
        self.max_column_dominance_spacer = max_column_dominance_spacer
        self.max_allowed_consecutive_spacers = max_allowed_consecutive_spacers
        self.max_allowed_same_spacers = max_allowed_same_spacers
        self.max_inconsistent_columns = max_inconsistent_columns
        self.min_avg_repeat_length = min_avg_repeat_length
        self.max_avg_repeat_length = max_avg_repeat_length
        self.max_avg_spacer_length = max_avg_spacer_length
        self.min_number_repeats = min_repeats

    @printing_if_filtered
    @exception_handler
    def _filter_by_column(self, candidate):
        def find_first_three_columns():
            list_three_columns = []
            list_gaped_repeats = candidate.list_gaped_repeats
            for index in range(len(list_gaped_repeats[0])):
                column_vec = [repeat[index] for repeat in list_gaped_repeats]
                column_gaps = sum([1 for x in column_vec if (x == " ")])
                percentage_gaps = column_gaps / len(column_vec)
                if percentage_gaps < 0.5:
                    list_three_columns.append(column_vec)
                if len(list_three_columns) == 3:
                    return list_three_columns

        def find_last_three_columns():
            list_three_columns = []
            list_gaped_repeats = candidate.list_gaped_repeats
            for index in range(len(list_gaped_repeats[0])-1, 0, -1):
                column_vec = [repeat[index] for repeat in list_gaped_repeats]
                column_gaps = sum([1 for x in column_vec if (x == " ")])
                percentage_gaps = column_gaps/len(column_vec)
                if percentage_gaps < 0.5:
                    list_three_columns.append(column_vec)
                if len(list_three_columns) == 3:
                    return list_three_columns

        for column in find_first_three_columns():
            column_characters = [x for x in column if (x not in (" ", "-"))]
            if column_characters:
                most_freq_char = max(column_characters, key=column_characters.count)
                most_freq_char_freq = column_characters.count(most_freq_char)
                freq = most_freq_char_freq/len(column_characters)
                if len(column) <= 4:
                    if freq < 0.49:
                        return False
                else:
                    if freq < self.column_dominance:
                        return False
            else:
                return False

        for column in find_last_three_columns():
            column_characters = [x for x in column if (x not in (" ", "-"))]
            if column_characters:
                most_freq_char = max(column_characters, key=column_characters.count)
                most_freq_char_freq = column_characters.count(most_freq_char)
                freq = most_freq_char_freq/len(column_characters)
                if len(column) <= 4:
                    if freq < 0.49:
                        return False
                else:
                    if freq < self.column_dominance:
                        return False
            else:
                return False
        return True

    @printing_if_filtered
    @exception_handler
    def _filter_by_min_avg_spacer(self, candidate):
        list_spacers = candidate.list_spacers
        avg_len = sum(len(x) for x in list_spacers) / len(list_spacers)
        if avg_len > self.min_avg_spacer_length:
            return True
        return False

    @printing_if_filtered
    @exception_handler
    def _filter_by_max_spacer(self, candidate):
        list_spacers = candidate.list_spacers
        long_spacers = [spacer for spacer in list_spacers if len(spacer) > self.max_spacer_length]
        if len(long_spacers) / len(list_spacers) > 0.3:
            return False
        if len(long_spacers) > 3:
            return False
        return True

    @printing_if_filtered
    @exception_handler
    def _filter_by_spacer_begin_end_similarity(self, candidate):
        list_spacers = candidate.list_spacers
        if len(list_spacers) >= 2:
            column_begin = [spacer[0] for spacer in list_spacers if spacer]
            most_freq_char_begin = max(column_begin, key=column_begin.count)
            most_freq_char_freq_begin = column_begin.count(most_freq_char_begin)

            freq_begin = most_freq_char_freq_begin / len(column_begin)
            if freq_begin > self.max_column_dominance_spacer:
                return False

            column_end = [spacer[-1] for spacer in list_spacers if spacer]
            most_freq_char_end = max(column_end, key=column_end.count)
            most_freq_char_freq_end = column_end.count(most_freq_char_end)

            freq_end = most_freq_char_freq_end / len(column_end)
            if freq_end > self.max_column_dominance_spacer:
                return False
        return True

    @printing_if_filtered
    @exception_handler
    def _filter_by_the_same_spacer(self, candidate):
        list_spacers = candidate.list_spacers
        list_spacers = [s for s in list_spacers if s]
        groups = [len(list(group)) for key, group in groupby(list_spacers)]
        if self.max_allowed_consecutive_spacers:
            if max(groups) > self.max_allowed_consecutive_spacers:
                return False

        list_sorted_spacers = sorted(list_spacers)
        groups_sorted = [len(list(group)) for key, group in groupby(list_sorted_spacers)]
        if self.max_allowed_same_spacers:
            if max(groups_sorted) > self.max_allowed_same_spacers:
                return False
        return True

    @printing_if_filtered
    @exception_handler
    def _filter_by_overall_repeat_consistency(self, candidate):
        list_column_consistency = []
        list_repeats_gaped = candidate.list_gaped_repeats
        for index, _ in enumerate(list_repeats_gaped[0]):
            column = [repeat[index] for repeat in list_repeats_gaped]
            column_characters = [x for x in column if (x not in (" ", "-"))]
            try:
                most_freq_char = max(column_characters, key=column_characters.count)
                most_freq_char_freq = column_characters.count(most_freq_char)
                freq = most_freq_char_freq / len(column_characters)
                list_column_consistency.append(freq)
            except ValueError:
                pass

        number_inconsistent = sum(1 for x in list_column_consistency if x < 0.66)
        if number_inconsistent > self.max_inconsistent_columns:
            return False
        return True

    @printing_if_filtered
    @exception_handler
    def _filter_min_number_repeats(self, candidate):
        list_repeats = candidate.list_repeats
        if len(list_repeats) >= self.min_number_repeats:
            return True
        return False

    @printing_if_filtered
    @exception_handler
    def _filter_min_avg_repeat_length(self, candidate):
        list_repeats = candidate.list_repeats
        avg_len = sum(len(x) for x in list_repeats) / len(list_repeats)
        if avg_len >= self.min_avg_repeat_length:
            return True
        return False

    @printing_if_filtered
    @exception_handler
    def _filter_max_avg_repeat_length(self, candidate):
        list_repeats = candidate.list_repeats
        avg_len = sum(len(x) for x in list_repeats) / len(list_repeats)
        if avg_len <= self.max_avg_repeat_length:
            return True
        return False

    @printing_if_filtered
    @exception_handler
    def _filter_max_avg_spacer_length(self, candidate):
        list_spacers = candidate.list_spacers
        if len(list_spacers) > 4:
            avg_len = sum(len(x) for x in list_spacers[1:-1]) / len(list_spacers)
            if avg_len <= self.max_avg_repeat_length:
                return True
        else:
            avg_len = sum(len(x) for x in list_spacers) / len(list_spacers)
            if avg_len <= self.max_avg_repeat_length:
                return True
        return False

    @printing_if_filtered
    @exception_handler
    def _filter_min_repeat_length(self, candidate):
        list_spacers = candidate.list_spacers
        avg_len = sum(len(x) for x in list_spacers) / len(list_spacers)
        if avg_len >= self.min_avg_repeat_length:
            return True
        return False

    def __call__(self, candidate):
        if not self._filter_by_column(candidate):
            return
        if not self._filter_by_min_avg_spacer(candidate):
            return
        if not self._filter_by_max_spacer(candidate):
            return
        if not self._filter_by_spacer_begin_end_similarity(candidate):
            return
        if not self._filter_by_the_same_spacer(candidate):
            return
        if not self._filter_by_overall_repeat_consistency(candidate):
            return
        if not self._filter_max_avg_repeat_length(candidate):
            return
        if not self._filter_min_avg_repeat_length(candidate):
            return
        if not self._filter_max_avg_spacer_length(candidate):
            return
        if not self._filter_min_number_repeats(candidate):
            return
        return candidate

#              CRISPR Candidate
#####################################################
#####################################################
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
            for rep_gapped in self.list_repeats_gaped:
                print(rep_gapped)
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

        self.list_gaped_repeats = self.list_repeats_gaped

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

    def dot_repr_web_server(self):
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
                spacer = "$" + self.list_spacers[index] + "$"
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

        string += "_" * 100 + "\n"

        string += "consensus: " + self.consensus + "\n"

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

    def __repr__(self):
        return self.dot_repr()

    def __eq__(self, other):
        if self.list_repeats == other.list_repeats:
            if self.list_repeats_gaped == other.list_repeats_gaped:
                if self.list_spacers == other.list_spacers:
                    return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)
