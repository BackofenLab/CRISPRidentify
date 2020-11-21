import os
from components_detection_refinement import CrisprCandidate


def rev_compliment_seq(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', ' ': ' ', "-": "-", ".": "."}
    try:
        compliment_seq = "".join([complement[nt] for nt in seq])
    except KeyError:
        compliment_seq = ""
        for nt in seq:
            if nt in complement:
                compliment_seq += complement[nt]
            else:
                compliment_seq += nt
    return compliment_seq[::-1]


class RevComComputation:
    def __init__(self, crispr_candidate):
        self.forward_crispr_candidate = crispr_candidate
        self.rev_com_candidate = None
        self._compute_rev_com_candidate()

    def _compute_rev_com_candidate(self):
        list_repeats_f = self.forward_crispr_candidate.list_repeats
        list_repeats_gaped_f = self.forward_crispr_candidate.list_repeats_gaped
        list_repeat_starts = self.forward_crispr_candidate.list_repeat_starts
        list_spacers = self.forward_crispr_candidate.list_spacers

        list_repeats_rev_com = [rev_compliment_seq(repeat) for repeat in list_repeats_f][::-1]

        list_repeats_gaped_rev_com = [rev_compliment_seq(repeat_gaped)
                                      for repeat_gaped in list_repeats_gaped_f][::-1]

        list_repeat_starts_rev_com = list_repeat_starts[::-1]

        list_spacers_rev_com = [spacer for spacer in list_spacers][::-1]

        self.rev_com_candidate = CrisprCandidate(list_repeats=list_repeats_rev_com,
                                                 list_repeats_gaped=list_repeats_gaped_rev_com,
                                                 list_repeat_starts=list_repeat_starts_rev_com,
                                                 list_spacers=list_spacers_rev_com)

    def output(self):
        return self.rev_com_candidate


class OutputCrispr:
    def __init__(self, crispr_candidate):
        stats = crispr_candidate.compute_stats()
        self.start = stats["start"]
        self.end = stats["end"]
        self.number_of_repeats = stats["number_repeats"]
        self.avg_repeat_length = stats["avg_repeat"]
        self.avg_spacer_length = stats["avg_spacer"]
        self.consensus = crispr_candidate.consensus
        self.list_repeats = crispr_candidate.list_repeats
        self.list_spacers = crispr_candidate.list_spacers

        self.list_repeat_starts = [x + 1 for x in crispr_candidate.list_repeat_starts]
        self.list_spacer_starts = [repeat_start + len(repeat) for repeat_start, repeat in zip(self.list_repeat_starts,
                                                                                              self.list_repeats)]
        self.dot_representation = crispr_candidate.dot_repr()
        self.dot_representation_web_server = crispr_candidate.dot_repr_web_server()


def index_crispr_candidates(self):
    for index, key in enumerate(sorted(self.dict_best.keys()), 1):
        self.dict_crispr_indexes[key] = index

    index = len(self.dict_best)
    cur_index = index + 1

    for key in sorted(self.dict_possible.keys()):
        if key not in self.dict_crispr_indexes:
            self.dict_crispr_indexes[key] = cur_index
            cur_index += 1

    for key in sorted(self.dict_possible_discarded.keys()):
        if key not in self.dict_crispr_indexes:
            self.dict_crispr_indexes[key] = cur_index
            cur_index += 1

    for key in sorted(self.dict_bad.keys()):
        if key not in self.dict_crispr_indexes:
            self.dict_crispr_indexes[key] = cur_index
            cur_index += 1


class SimpleOutputMaker:
    def __init__(self, categories, result_path, list_features):
        self.result_path = result_path
        self.categories = categories
        self.list_features = list_features

        self.dict_crispr_indexes = {}

        self._index_crispr_candidates()
        self._write_simple_txt_files()

    def _index_crispr_candidates(self):
        for index, key in enumerate(sorted(self.categories[0].keys()), 1):
            self.dict_crispr_indexes[key] = index

        index = len(self.categories[0])
        cur_index = index + 1

        for key in sorted(self.categories[2].keys()):
            if key not in self.dict_crispr_indexes:
                self.dict_crispr_indexes[key] = cur_index
                cur_index += 1

        for key in sorted(self.categories[3].keys()):
            if key not in self.dict_crispr_indexes:
                self.dict_crispr_indexes[key] = cur_index
                cur_index += 1

        for key in sorted(self.categories[4].keys()):
            if key not in self.dict_crispr_indexes:
                self.dict_crispr_indexes[key] = cur_index
                cur_index += 1
        
    def _write_simple_txt_files(self):
        if not os.path.exists(self.result_path):
            os.makedirs(self.result_path)

        f_bona_fide = open(self.result_path + '/Bona_Fide_Candidates.txt', 'w')
        f_alternative = open(self.result_path + '/Alternative_Candidates.txt', 'w')
        f_possible = open(self.result_path + '/Possible_Candidates.txt', 'w')
        f_possible_discarded = open(self.result_path + '/Possible_Discarded_Candidates.txt', 'w')
        f_low_score = open(self.result_path + '/Low_Score_Candidates.txt', 'w')

        file_names = [f_bona_fide, f_alternative, f_possible, f_possible_discarded, f_low_score]
        category_names = ["", "Alternative ", "Possible ", "Possible Discarded ", "Low score "]

        for index, file_name, category_name in zip(range(5), file_names, category_names):
            dictionary = self.categories[index]
            for key in sorted(dictionary):
                for candidate in dictionary[key]:
                    score = candidate[0]
                    crispr = candidate[1]
                    crispr_stats = crispr.compute_stats()
                    file_name.write("{}CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                                 .format(category_name, self.dict_crispr_indexes[key], crispr_stats["start"],
                                         crispr_stats["end"], crispr_stats["number_repeats"],
                                         crispr_stats["avg_repeat"], crispr_stats["avg_spacer"]))
    
                    file_name.write(crispr.dot_repr())
    
                    file_name.write("\n")
                    list_reported_features = []
                    for feature_index, feature_list in enumerate(self.list_features):
                        for feature, value in zip(feature_list, candidate[2][feature_index][0]):
                            if feature not in list_reported_features:
                                file_name.write("{}: {}\n".format(feature, value))
                                list_reported_features.append(feature)
    
                    file_name.write("\n")
                    file_name.write("Certainty Score: {}\n\n".format(score))
                file_name.write('\n{}\n\n'.format('=' * 100))
            file_name.close()


class SummaryOutputMaker:
    def __init__(self, result_path, categories, non_array_data, header, list_feature_names):
        self.result_path = result_path
        self.categories = categories
        self.non_array_data = non_array_data
        self.header = header
        self.list_feature_names = list_feature_names

        self._make_text_summary()

    def _make_text_summary(self):
        result_path = self.result_path + '/Summary.txt'
        list_crisprs = [list_info[0][1] for list_info in self.categories[0].values()]
        list_scores = [list_info[0][0] for list_info in self.categories[0].values()]
        list_feature_vectors = [list_info[0][2] for list_info in self.categories[0].values()]
        
        with open(result_path, "w") as f:
            f.write(self.header)
            f.write("\n")

            for index, array in enumerate(list_crisprs):
                if index in self.non_array_data["Cas"]:
                    cas_genes = self.non_array_data["Cas"][index]
                    if cas_genes:
                        f.write("Cas genes: ")
                        for index_cas, cluster in enumerate(cas_genes):
                            if index_cas != (len(cas_genes) - 1):
                                f.write("{} [{}-{}], ".format(cluster[2], cluster[0], cluster[1]))
                            else:
                                f.write("{} [{}-{}]".format(cluster[2], cluster[0], cluster[1]))
                    f.write("\n\n")

                strand = self.non_array_data["Strand"][index]
                if strand == "Forward":
                    output_crispr = OutputCrispr(array)
                else:
                    output_crispr = OutputCrispr(RevComComputation(array).output())

                crispr_index = index + 1
                start, end = output_crispr.start, output_crispr.end
                number_of_repeats = output_crispr.number_of_repeats
                avg_length_repeat, avg_length_spacer = output_crispr.avg_repeat_length, output_crispr.avg_spacer_length

                if index > 0:
                    f.write('\n{}\n\n'.format('=' * 100))
                f.write(
                    "CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                    .format(str(crispr_index), start, end, number_of_repeats,
                            avg_length_repeat, avg_length_spacer))

                f.write(output_crispr.dot_representation_web_server)

                f.write("\n")

                f.write("Leader region\n")

                if strand == "Reverse":
                    leader = self.non_array_data["Leader"][index]
                else:
                    leader = self.non_array_data["Leader_rev_com"][index]

                f.write(leader)

                f.write("\n\nStrand: {}\n\n".format(strand))

                f.write("#   Array features:\n")
                list_reported_features = []
                for group_of_features, feature_vector in zip(self.list_feature_names,
                                                             list_feature_vectors[index]):
                    for feature_name, feature_value in zip(group_of_features, feature_vector[0]):
                        if feature_name not in list_reported_features:
                            f.write("#   {}: {}\n".format(feature_name, feature_value))
                            list_reported_features.append(feature_name)

                f.write("_" * 30)
                f.write("\n")
                f.write("#   Certainty Score: {}\n\n".format(list_scores[index]))

                if index in self.non_array_data["IS"]:
                    f.write("IS Element: {} [{}-{}]\n\n".format(self.non_array_data["IS"][index][4],
                                                                self.non_array_data["IS"][index][0],
                                                                self.non_array_data["IS"][index][1]))

            last_index = len(list_crisprs)
            if last_index in self.non_array_data["Cas"]:
                f.write("\n\n")
                cas_genes = self.non_array_data["Cas"][last_index]
                if cas_genes:
                    f.write("Cas genes: ")
                    for index, cluster in enumerate(cas_genes):
                        if index != (len(cas_genes) - 1):
                            f.write("{} [{}-{}], ".format(cluster[2], cluster[0], cluster[1]))
                        else:
                            f.write("{} [{}-{}]".format(cluster[2], cluster[0], cluster[1]))
                f.write("\n\n")