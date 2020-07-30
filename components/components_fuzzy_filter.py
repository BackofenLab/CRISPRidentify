from functools import wraps


DEBUG_MODE = True


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
                #print("\n\n")
                #print(args)
                #print(function.__name__)
                #print("\n\n")

                with open("filtered_results.txt", "a") as f:
                    f.write("\n\n")
                    f.write(str(args))
                    f.write("\n\n")
                    f.write(function.__name__)
                    f.write("\n\n")

        return result
    return wrapper


class AdvancedFuzzySearchFilter:
    def __init__(self, min_column_dominance_repeat, min_avg_spacer_length,
                 max_spacer_length, max_column_dominance_spacer, max_allowed_conseq_spacers,
                 max_allowed_same_spacers, max_inconsistent_columns, min_avg_repeat_length,
                 max_avg_repeat_length, max_avg_spacer_length, min_repeats):

        self.column_dominance = min_column_dominance_repeat
        self.min_avg_spacer_length = min_avg_spacer_length
        self.max_spacer_length = max_spacer_length
        self.max_column_dominance_spacer = max_column_dominance_spacer
        self.max_allowed_conseq_spacers = max_allowed_conseq_spacers
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
            most_freq_char = max(column_characters, key=column_characters.count)
            most_freq_char_freq = column_characters.count(most_freq_char)
            freq = most_freq_char_freq/len(column_characters)
            if len(column) == 4:
                if freq < 0.49:
                    return False
            else:
                if freq < self.column_dominance:
                    return False

        for column in find_last_three_columns():
            column_characters = [x for x in column if (x not in (" ", "-"))]
            most_freq_char = max(column_characters, key=column_characters.count)
            most_freq_char_freq = column_characters.count(most_freq_char)
            freq = most_freq_char_freq/len(column_characters)
            if len(column) == 4:
                if freq < 0.49:
                    return False
            else:
                if freq < self.column_dominance:
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
        if len(list_spacers) > 2:
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
        """
        list_spacers = candidate.list_spacers
        groups = [len(list(group)) for key, group in groupby(list_spacers)]
        if max(groups) > self.max_allowed_conseq_spacers:
            return False

        list_sorted_spacers = sorted(list_spacers)
        groups_sorted = [len(list(group)) for key, group in groupby(list_sorted_spacers)]
        if max(groups_sorted) > self.max_allowed_same_spacers:
            return False
        """
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
        if not self._filter_by_min_avg_spacer(candidate):
            return
        if not self._filter_min_number_repeats(candidate):
            return
        return candidate