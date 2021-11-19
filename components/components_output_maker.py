import os
import pickle
import json
from os import listdir
from os.path import basename
from os.path import isfile, join
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

        list_spacers_rev_com = [rev_compliment_seq(spacer) for spacer in list_spacers][::-1]

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


class SimpleOutputMaker:
    def __init__(self, categories, non_array_data, result_path, list_features):
        self.result_path = result_path
        self.categories = categories
        self.non_array_data = non_array_data
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

        indexes_bona_fide = [self.dict_crispr_indexes[key] for key in sorted(self.categories[0].keys())]
        indexes_alternative = [self.dict_crispr_indexes[key] for key in sorted(self.categories[1].keys())]
        indexes_possible = [self.dict_crispr_indexes[key] for key in sorted(self.categories[2].keys())]
        indexes_possible_discarded = [self.dict_crispr_indexes[key] for key in sorted(self.categories[3].keys())]
        indexes_low_score = [self.dict_crispr_indexes[key] for key in sorted(self.categories[4].keys())]

        self.list_indexes = [indexes_bona_fide, indexes_alternative, indexes_possible,
                             indexes_possible_discarded, indexes_low_score]

    def _write_simple_txt_files(self):
        if not os.path.exists(self.result_path):
            os.makedirs(self.result_path)

        f_alternative = open(self.result_path + '/Alternative_Candidates.txt', 'w')
        f_possible = open(self.result_path + '/Possible_Candidates.txt', 'w')
        f_possible_discarded = open(self.result_path + '/Possible_Discarded_Candidates.txt', 'w')
        f_low_score = open(self.result_path + '/Low_Score_Candidates.txt', 'w')

        file_names = [f_alternative, f_possible, f_possible_discarded, f_low_score]
        category_names = ["Alternative", "Possible", "Possible Discarded", "Low score"]

        for category_index, category_name, file_name in zip(range(1, 5),
                                                            category_names,
                                                            file_names):

            arrays = [el[1] for key in self.categories[category_index].keys()
                      for el in self.categories[category_index][key]]
            scores = [el[0] for key in self.categories[category_index].keys()
                      for el in self.categories[category_index][key]]

            features = [el[2] for key in self.categories[category_index].keys()
                        for el in self.categories[category_index][key]]

            array_indexes = self.list_indexes[category_index]

            if array_indexes:
                for index, array_index, array, score, feature_info in zip(range(len(arrays)), array_indexes,
                                                                          arrays, scores, features):
                    if category_name in self.non_array_data["Strand"]:
                        strand = self.non_array_data["Strand"][category_name][index]
                    else:
                        strand = "Forward (Orientation was not computed)"
                    if strand == "Reversed":
                        crispr = RevComComputation(array).output()
                    else:
                        crispr = array

                    crispr_stats = crispr.compute_stats()
                    file_name.write(
                        "{} CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                            .format(category_name, array_index, crispr_stats["start"],
                                    crispr_stats["end"], crispr_stats["number_repeats"],
                                    crispr_stats["avg_repeat"], crispr_stats["avg_spacer"]))

                    file_name.write(crispr.dot_repr())

                    file_name.write(f"\nStrand: {strand}\n\n")

                    file_name.write("\n")
                    list_reported_features = []
                    for feature_index, feature_list in enumerate(self.list_features):
                        for feature, value in zip(feature_list, feature_info[feature_index][0]):
                            if feature not in list_reported_features:
                                file_name.write("{}: {}\n".format(feature, value))
                                list_reported_features.append(feature)

                    file_name.write("\n")
                    file_name.write("Certainty Score: {}\n\n\n".format(score))

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
        result_path = self.result_path + '/Bona-Fide_Candidates.txt'
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

                strand = self.non_array_data["Strand"]["Bona-fide"][index]
                if strand in ("Forward", "Forward (Orientation was not computed)"):
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

                f.write(output_crispr.dot_representation)

                f.write("\n")

                f.write("Leader region\n")

                leader = self.non_array_data["Leader"][0][index]

                f.write(leader)

                f.write("\n\n")

                f.write("Downstream region\n")

                downstream = self.non_array_data["Downstream"][0][index]

                f.write(downstream)

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


class SummaryMakerCSV:
    def __init__(self, result_path, categories, non_array_data):
        self.result_path = result_path
        self.categories = categories
        self.non_array_data = non_array_data

        self.dict_crispr_indexes = {}

        self._index_crispr_candidates()
        self._write_csv_summary()

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

        indexes_bona_fide = [self.dict_crispr_indexes[key] for key in sorted(self.categories[0].keys())]
        indexes_alternative = [self.dict_crispr_indexes[key] for key in sorted(self.categories[1].keys())]
        indexes_possible = [self.dict_crispr_indexes[key] for key in sorted(self.categories[2].keys())]
        self.list_indexes = [indexes_bona_fide, indexes_alternative, indexes_possible]

    def _write_csv_summary(self):

        all_arrays = [el[1] for category_index in range(3) for key in self.categories[category_index].keys()
                      for el in self.categories[category_index][key]]

        if all_arrays:

            result_csv_path = self.result_path + '/Summary.csv'
            with open(result_csv_path, "w") as f:
                f.write(",".join(["ID", "Region index", "Start", "End", "Length", "Consensus repeat", "Repeat Length",
                                  "Average Spacer Length", "Number of spacers", "Strand", "Category"]))
                f.write("\n")
                global_index = 1
                for category_index, category in zip(range(3), ["Bona-fide", "Alternative", "Possible"]):
                    arrays = [el[1] for key in self.categories[category_index].keys()
                              for el in self.categories[category_index][key]]
                    array_indexes = self.list_indexes[category_index]
                    for index, array_index, array in zip(range(len(arrays)), array_indexes, arrays):
                        strand = self.non_array_data["Strand"][category][index]
                        crispr = array
                        crispr_stats = crispr.compute_stats()
                        crispr_index = str(array_index)
                        start = str(crispr_stats["start"])
                        end = str(crispr_stats["end"])
                        length = str(int(end) - int(start) + 1)
                        if strand in ("Forward", "Forward (Orientation was not computed)"):
                            consensus_repeat = crispr.consensus
                        else:
                            consensus_repeat = rev_compliment_seq(crispr.consensus)

                        repeat_length = str(crispr_stats["avg_repeat"])
                        average_spacer_length = str(crispr_stats["avg_spacer"])
                        number_of_spacers = str(crispr_stats["number_repeats"] - 1)

                        string_to_write = ",".join([str(global_index), crispr_index, start, end, length,
                                                    consensus_repeat, repeat_length, average_spacer_length,
                                                    number_of_spacers, strand, category])
                        global_index += 1
                        f.write(string_to_write)
                        f.write("\n")


class PickleOutputMaker:
    def __init__(self, file_path, pickle_result_folder, parameters, categories,
                 non_array_data, header, list_feature_names):
        self.file_name = file_path
        self.pickle_result_folder = pickle_result_folder
        self.parameters = parameters
        self.categories = categories
        self.non_array_data = non_array_data
        self.header = header
        self.list_feature_names = list_feature_names

        self._write_pickle()

    def _write_pickle(self):
        try:
            os.mkdir(self.pickle_result_folder)
        except OSError:
            pass

        file_base = basename(self.file_name)
        acc_num = file_base.split(".")[0]
        pickle.dump(self.categories, open(self.pickle_result_folder + '/' + acc_num + '.pkl', "wb"))


class JsonOutputMaker:
    def __init__(self, file_path, json_result_folder, categories, non_array_data, list_feature_names):
        self.file_path = file_path
        self.json_result_folder = json_result_folder
        self.categories = categories
        self.non_array_data = non_array_data
        self.list_feature_names = list_feature_names

        self.dict_crispr_indexes = {}

        self._index_crispr_candidates()
        self._write_json()

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

        indexes_bona_fide = [self.dict_crispr_indexes[key] for key in sorted(self.categories[0].keys())]
        indexes_alternative = [self.dict_crispr_indexes[key] for key in sorted(self.categories[1].keys())]
        indexes_possible = [self.dict_crispr_indexes[key] for key in sorted(self.categories[2].keys())]
        indexes_possible_discarded = [self.dict_crispr_indexes[key] for key in sorted(self.categories[3].keys())]
        indexes_low_score = [self.dict_crispr_indexes[key] for key in sorted(self.categories[4].keys())]
        self.list_indexes = [indexes_bona_fide, indexes_alternative, indexes_possible,
                             indexes_possible_discarded, indexes_low_score]

    def _write_json(self):
        try:
            os.mkdir(self.json_result_folder)
        except OSError:
            pass
        json_dictionary = {}
        global_index = 1
        dictionary_all_arrays = {}

        category_names = ["Bona-fide", "Alternative", "Possible", "Possible Discarded", "Low score"]

        for category_index, category_name in zip(range(0, 5),
                                                 category_names):
            dict_category = {}

            arrays = [el[1] for key in self.categories[category_index].keys()
                      for el in self.categories[category_index][key]]
            scores = [el[0] for key in self.categories[category_index].keys()
                      for el in self.categories[category_index][key]]

            features = [el[2] for key in self.categories[category_index].keys()
                        for el in self.categories[category_index][key]]

            array_indexes = self.list_indexes[category_index]

            if array_indexes:
                for index, array_index, array, score, feature_info in zip(range(len(arrays)), array_indexes,
                                                                          arrays, scores, features):
                    if category_name in self.non_array_data["Strand"]:
                        strand = self.non_array_data["Strand"][category_name][index]
                    else:
                        strand = "Forward (Orientation was not computed)"
                    if strand == "Reversed":
                        crispr = RevComComputation(array).output()
                    else:
                        crispr = array

                    crispr_stats = crispr.compute_stats()

                    crispr_start = crispr_stats["start"]
                    crispr_end = crispr_stats["end"] + 1
                    crispr_length = crispr_end - crispr_start + 1

                    dict_crispr = self.crispr_candidate_to_dictionary(crispr)
                    dict_crispr["strand"] = strand
                    dict_crispr["start"] = crispr_start
                    dict_crispr["end"] = crispr_end
                    dict_crispr["length"] = crispr_length
                    dict_crispr["array_index"] = array_index
                    dict_crispr["certainty_score"] = score

                    dict_category[global_index] = dict_crispr
                    global_index += 1

            dictionary_all_arrays[category_name] = dict_category

        file_base = basename(self.file_path)
        acc_num = file_base.split(".")[0]

        json_dictionary["Arrays"] = dictionary_all_arrays

        json_representation = json.dumps(json_dictionary)
        with open(self.json_result_folder + '/' + acc_num + '.json', "w") as f:
            f.write(json_representation)

    @staticmethod
    def crispr_candidate_to_dictionary(crispr_candidate):
        list_repeat_starts = crispr_candidate.list_repeat_starts
        list_repeats = crispr_candidate.list_repeats
        list_spacers = crispr_candidate.list_spacers
        consensus_repeat = crispr_candidate.consensus
        dot_representation = crispr_candidate.dot_repr()
        dot_representation_web_server = crispr_candidate.dot_repr_web_server()

        dict_crispr = {"list_repeats": list_repeats,
                       "list_spacers": list_spacers,
                       "list_repeat_starts": list_repeat_starts,
                       "consensus_repeat": consensus_repeat,
                       "dot_representation": dot_representation,
                       "dot_representation_web_server": dot_representation_web_server}

        return dict_crispr


class GFFOutputMaker:
    def __init__(self, result_path, categories, non_array_data, header, list_feature_names):
        self.result_path = result_path
        self.categories = categories
        self.non_array_data = non_array_data
        self.list_features = list_feature_names

        self.acc_num = result_path.split("/")[-1]
        if not self.acc_num:
            self.acc_num = result_path.split("/")[-2]

        self.dict_crispr_indexes = {}

        self._index_crispr_candidates()
        self._create_gff_files()

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

        indexes_bona_fide = [self.dict_crispr_indexes[key] for key in sorted(self.categories[0].keys())]
        indexes_alternative = [self.dict_crispr_indexes[key] for key in sorted(self.categories[1].keys())]
        indexes_possible = [self.dict_crispr_indexes[key] for key in sorted(self.categories[2].keys())]
        self.list_indexes = [indexes_bona_fide, indexes_alternative, indexes_possible]

    def _create_gff_files(self):
        self.gff_folder = join(self.result_path, "gff_result")
        directory_exists = os.path.isdir(self.gff_folder)
        if not directory_exists:
            command = "mkdir " + self.gff_folder
            try:
                os.system(command)
            except Exception:
                pass

        self._create_gff_bona_fide()
        self._create_gff_alternative()
        self._create_gff_possible()
        self._create_gff_complete()

    def _create_gff_bona_fide(self):
        with open(join(self.gff_folder, "bona-fide.gff"), "w") as f:
            arrays = [el[1] for key in self.categories[0].keys()
                      for el in self.categories[0][key]]
            scores = [el[0] for key in self.categories[0].keys()
                      for el in self.categories[0][key]]

            features = [el[2] for key in self.categories[0].keys()
                        for el in self.categories[0][key]]

            array_indexes = self.list_indexes[0]

            if array_indexes:
                for index, array_index, array, score, feature_info in zip(range(len(arrays)), array_indexes,
                                                                          arrays, scores, features):
                    if "Bona-fide" in self.non_array_data["Strand"]:
                        strand = self.non_array_data["Strand"]["Bona_fide"][index]
                    else:
                        strand = "Forward (Orientation was not computed)"

                    strand_sign = "+" if strand in ["Forward", "Forward (Orientation was not computed)"] else "-"

                    crispr = array

                    crispr_stats = crispr.compute_stats()
                    crispr_start = crispr_stats["start"]
                    crispr_end = crispr_stats["end"] + 1
                    crispr_length = crispr_end - crispr_start + 1
                    consensus = crispr.consensus
                    line_crispr = f"{self.acc_num}\tCRISPRidentify\tbona-fide_array_region\t{crispr_start}\t{crispr_end}\t{crispr_length}\t{strand_sign}\t.\tID=CRISPR{array_index}_{crispr_start}_{crispr_end};Note={consensus};Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score={score}\n"
                    f.write(line_crispr)

                    repeats = crispr.list_repeats
                    spacers = crispr.list_spacers
                    repeat_starts = crispr.list_repeat_starts

                    repeat_starts = [1 + r_s for r_s in repeat_starts]
                    repeat_ends = [rs + len(repeat) for rs, repeat in zip(repeat_starts, repeats)]

                    spacers_starts = [r_e + 1 for r_e in repeat_ends[:-1]]
                    spacers_ends = [s_s + len(s) for s_s, s in zip(spacers_starts, spacers)]

                    repeat_indexes = list(range(1, len(repeats) + 1))
                    spacer_indexes = list(range(1, len(spacers) + 1))

                    for r, s, r_s, r_e, s_s, s_e, r_index, s_index in zip(repeats, spacers, repeat_starts,
                                                                          repeat_ends, spacers_starts, spacers_ends,
                                                                          repeat_indexes, spacer_indexes):
                        repeat_len = len(r)
                        spacer_len = len(s)
                        line_repeat = f"{self.acc_num}\tCRISPRidentify\tdirect_repeat\t{r_s}\t{r_e}\t{repeat_len}\t{strand_sign}\t.\tID=CRISPR{array_index}_REPEAT{r_index}_{r_s}_{r_e};Name=CRISPR{array_index}_REPEAT{r_index}_{r_s}_{r_e};Parent=CRISPR{array_index}_{r_s}_{r_e};Note={r};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                        line_spacer = f"{self.acc_num}\tCRISPRidentify\tspacer\t{s_s}\t{s_e}\t{spacer_len}\t{strand_sign}\t.\tID=CRISPR{array_index}_SPACER{s_index}_{s_s}_{s_e};Name=CRISPR{array_index}_REPEAT{s_index}_{s_s}_{s_e};Parent=CRISPR{array_index}_{s_s}_{s_e};Note={r};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                        f.write(line_repeat)
                        f.write(line_spacer)

                    last_repeat = repeats[-1]
                    len_last_repeat = len(last_repeat)
                    last_repeat_start = repeat_starts[-1]
                    last_repeat_end = repeat_ends[-1]
                    last_repeat_index = repeat_indexes[-1]

                    line_repeat = f"{self.acc_num}\tCRISPRidentify\tdirect_repeat\t{last_repeat_start}\t{last_repeat_end}\t{len_last_repeat}\t{strand_sign}\t.\tID=CRISPR{array_index}_REPEAT{last_repeat_index}_{last_repeat_start}_{last_repeat_end};Name=CRISPR{array_index}_REPEAT{array_index}_{last_repeat_start}_{last_repeat_end};Parent=CRISPR{array_index}_{last_repeat_start}_{last_repeat_end};Note={last_repeat};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                    f.write(line_repeat)

    def _create_gff_alternative(self):
        with open(join(self.gff_folder, "alternative.gff"), "w") as f:
            arrays = [el[1] for key in self.categories[1].keys()
                      for el in self.categories[1][key]]
            scores = [el[0] for key in self.categories[1].keys()
                      for el in self.categories[1][key]]

            features = [el[2] for key in self.categories[1].keys()
                        for el in self.categories[1][key]]

            array_indexes = self.list_indexes[1]

            if array_indexes:
                for index, array_index, array, score, feature_info in zip(range(len(arrays)), array_indexes,
                                                                          arrays, scores, features):
                    if "Alternative" in self.non_array_data["Strand"]:
                        strand = self.non_array_data["Strand"]["Alternative"][index]
                    else:
                        strand = "Forward (Orientation was not computed)"

                    strand_sign = "+" if strand in ["Forward", "Forward (Orientation was not computed)"] else "-"

                    crispr = array

                    crispr_stats = crispr.compute_stats()
                    crispr_start = crispr_stats["start"]
                    crispr_end = crispr_stats["end"] + 1
                    crispr_length = crispr_end - crispr_start + 1
                    consensus = crispr.consensus
                    line_crispr = f"{self.acc_num}\tCRISPRidentify\talternative_array_region\t{crispr_start}\t{crispr_end}\t{crispr_length}\t{strand_sign}\t.\tID=CRISPR{array_index}_{crispr_start}_{crispr_end};Note={consensus};Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score={score}\n"
                    f.write(line_crispr)

                    repeats = crispr.list_repeats
                    spacers = crispr.list_spacers
                    repeat_starts = crispr.list_repeat_starts

                    repeat_starts = [1 + r_s for r_s in repeat_starts]
                    repeat_ends = [rs + len(repeat) for rs, repeat in zip(repeat_starts, repeats)]

                    spacers_starts = [r_e + 1 for r_e in repeat_ends[:-1]]
                    spacers_ends = [s_s + len(s) for s_s, s in zip(spacers_starts, spacers)]

                    repeat_indexes = list(range(1, len(repeats) + 1))
                    spacer_indexes = list(range(1, len(spacers) + 1))

                    for r, s, r_s, r_e, s_s, s_e, r_index, s_index in zip(repeats, spacers, repeat_starts,
                                                                          repeat_ends, spacers_starts, spacers_ends,
                                                                          repeat_indexes, spacer_indexes):
                        repeat_len = len(r)
                        spacer_len = len(s)
                        line_repeat = f"{self.acc_num}\tCRISPRidentify\tdirect_repeat\t{r_s}\t{r_e}\t{repeat_len}\t{strand_sign}\t.\tID=CRISPR{array_index}_REPEAT{r_index}_{r_s}_{r_e};Name=CRISPR{array_index}_REPEAT{r_index}_{r_s}_{r_e};Parent=CRISPR{array_index}_{r_s}_{r_e};Note={r};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                        line_spacer = f"{self.acc_num}\tCRISPRidentify\tspacer\t{s_s}\t{s_e}\t{spacer_len}\t{strand_sign}\t.\tID=CRISPR{array_index}_SPACER{s_index}_{s_s}_{s_e};Name=CRISPR{array_index}_REPEAT{s_index}_{s_s}_{s_e};Parent=CRISPR{array_index}_{s_s}_{s_e};Note={r};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                        f.write(line_repeat)
                        f.write(line_spacer)

                    last_repeat = repeats[-1]
                    len_last_repeat = len(last_repeat)
                    last_repeat_start = repeat_starts[-1]
                    last_repeat_end = repeat_ends[-1]
                    last_repeat_index = repeat_indexes[-1]

                    line_repeat = f"{self.acc_num}\tCRISPRidentify\tdirect_repeat\t{last_repeat_start}\t{last_repeat_end}\t{len_last_repeat}\t{strand_sign}\t.\tID=CRISPR{array_index}_REPEAT{last_repeat_index}_{last_repeat_start}_{last_repeat_end};Name=CRISPR{array_index}_REPEAT{array_index}_{last_repeat_start}_{last_repeat_end};Parent=CRISPR{array_index}_{last_repeat_start}_{last_repeat_end};Note={last_repeat};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                    f.write(line_repeat)

    def _create_gff_possible(self):
        with open(join(self.gff_folder, "possible.gff"), "w") as f:
            arrays = [el[1] for key in self.categories[2].keys()
                      for el in self.categories[2][key]]
            scores = [el[0] for key in self.categories[2].keys()
                      for el in self.categories[2][key]]

            features = [el[2] for key in self.categories[2].keys()
                        for el in self.categories[2][key]]

            array_indexes = self.list_indexes[2]

            if array_indexes:
                for index, array_index, array, score, feature_info in zip(range(len(arrays)), array_indexes,
                                                                          arrays, scores, features):
                    if "Possible" in self.non_array_data["Strand"]:
                        strand = self.non_array_data["Strand"]["Possible"][index]
                    else:
                        strand = "Forward (Orientation was not computed)"

                    strand_sign = "+" if strand in ["Forward", "Forward (Orientation was not computed)"] else "-"

                    crispr = array

                    crispr_stats = crispr.compute_stats()
                    crispr_start = crispr_stats["start"]
                    crispr_end = crispr_stats["end"] + 1
                    crispr_length = crispr_end - crispr_start + 1
                    consensus = crispr.consensus
                    line_crispr = f"{self.acc_num}\tCRISPRidentify\tpossible_array_region\t{crispr_start}\t{crispr_end}\t{crispr_length}\t{strand_sign}\t.\tID=CRISPR{array_index}_{crispr_start}_{crispr_end};Note={consensus};Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score={score}\n"
                    f.write(line_crispr)

                    repeats = crispr.list_repeats
                    spacers = crispr.list_spacers
                    repeat_starts = crispr.list_repeat_starts

                    repeat_starts = [1 + r_s for r_s in repeat_starts]
                    repeat_ends = [rs + len(repeat) for rs, repeat in zip(repeat_starts, repeats)]

                    spacers_starts = [r_e + 1 for r_e in repeat_ends[:-1]]
                    spacers_ends = [s_s + len(s) for s_s, s in zip(spacers_starts, spacers)]

                    repeat_indexes = list(range(1, len(repeats) + 1))
                    spacer_indexes = list(range(1, len(spacers) + 1))

                    for r, s, r_s, r_e, s_s, s_e, r_index, s_index in zip(repeats, spacers, repeat_starts,
                                                                          repeat_ends, spacers_starts, spacers_ends,
                                                                          repeat_indexes, spacer_indexes):
                        repeat_len = len(r)
                        spacer_len = len(s)
                        line_repeat = f"{self.acc_num}\tCRISPRidentify\tdirect_repeat\t{r_s}\t{r_e}\t{repeat_len}\t{strand_sign}\t.\tID=CRISPR{array_index}_REPEAT{r_index}_{r_s}_{r_e};Name=CRISPR{array_index}_REPEAT{r_index}_{r_s}_{r_e};Parent=CRISPR{array_index}_{r_s}_{r_e};Note={r};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                        line_spacer = f"{self.acc_num}\tCRISPRidentify\tspacer\t{s_s}\t{s_e}\t{spacer_len}\t{strand_sign}\t.\tID=CRISPR{array_index}_SPACER{s_index}_{s_s}_{s_e};Name=CRISPR{array_index}_REPEAT{s_index}_{s_s}_{s_e};Parent=CRISPR{array_index}_{s_s}_{s_e};Note={r};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                        f.write(line_repeat)
                        f.write(line_spacer)

                    last_repeat = repeats[-1]
                    len_last_repeat = len(last_repeat)
                    last_repeat_start = repeat_starts[-1]
                    last_repeat_end = repeat_ends[-1]
                    last_repeat_index = repeat_indexes[-1]

                    line_repeat = f"{self.acc_num}\tCRISPRidentify\tdirect_repeat\t{last_repeat_start}\t{last_repeat_end}\t{len_last_repeat}\t{strand_sign}\t.\tID=CRISPR{array_index}_REPEAT{last_repeat_index}_{last_repeat_start}_{last_repeat_end};Name=CRISPR{array_index}_REPEAT{array_index}_{last_repeat_start}_{last_repeat_end};Parent=CRISPR{array_index}_{last_repeat_start}_{last_repeat_end};Note={last_repeat};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                    f.write(line_repeat)

    def _create_gff_complete(self):
        with open(join(self.gff_folder, "combined.gff"), "w") as f:

            for category_index, category_name in zip([0, 1, 2], ["Bona-fide", "Alternative", "Possible"]):

                arrays = [el[1] for key in self.categories[category_index].keys()
                          for el in self.categories[category_index][key]]
                scores = [el[0] for key in self.categories[category_index].keys()
                          for el in self.categories[category_index][key]]

                features = [el[2] for key in self.categories[category_index].keys()
                            for el in self.categories[category_index][key]]

                array_indexes = self.list_indexes[category_index]

                if array_indexes:
                    for index, array_index, array, score, feature_info in zip(range(len(arrays)), array_indexes,
                                                                              arrays, scores, features):
                        if category_name in self.non_array_data["Strand"]:
                            strand = self.non_array_data["Strand"][category_name][index]
                        else:
                            strand = "Forward (Orientation was not computed)"

                        strand_sign = "+" if strand in ["Forward", "Forward (Orientation was not computed)"] else "-"

                        crispr = array

                        crispr_stats = crispr.compute_stats()
                        crispr_start = crispr_stats["start"]
                        crispr_end = crispr_stats["end"]
                        crispr_length = crispr_end - crispr_start + 1
                        consensus = crispr.consensus
                        line_crispr = f"{self.acc_num}\tCRISPRidentify\t{category_name}_array_region\t{crispr_start}\t{crispr_end}\t{crispr_length}\t{strand_sign}\t.\tID=CRISPR{array_index}_{crispr_start}_{crispr_end};Note={consensus};Dbxref=SO:0001459;Ontology_term=CRISPR;Array_quality_score={score}\n"
                        f.write(line_crispr)

                        repeats = crispr.list_repeats
                        spacers = crispr.list_spacers
                        repeat_starts = crispr.list_repeat_starts

                        repeat_starts = [1 + r_s for r_s in repeat_starts]
                        repeat_ends = [rs + len(repeat) for rs, repeat in zip(repeat_starts, repeats)]

                        spacers_starts = [r_e + 1 for r_e in repeat_ends[:-1]]
                        spacers_ends = [s_s + len(s) for s_s, s in zip(spacers_starts, spacers)]

                        repeat_indexes = list(range(1, len(repeats) + 1))
                        spacer_indexes = list(range(1, len(spacers) + 1))

                        for r, s, r_s, r_e, s_s, s_e, r_index, s_index in zip(repeats, spacers, repeat_starts,
                                                                              repeat_ends, spacers_starts, spacers_ends,
                                                                              repeat_indexes, spacer_indexes):
                            repeat_len = len(r)
                            spacer_len = len(s)
                            line_repeat = f"{self.acc_num}\tCRISPRidentify\tdirect_repeat\t{r_s}\t{r_e}\t{repeat_len}\t{strand_sign}\t.\tID=CRISPR{array_index}_REPEAT{r_index}_{r_s}_{r_e};Name=CRISPR{array_index}_REPEAT{r_index}_{r_s}_{r_e};Parent=CRISPR{array_index}_{r_s}_{r_e};Note={r};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                            line_spacer = f"{self.acc_num}\tCRISPRidentify\tspacer\t{s_s}\t{s_e}\t{spacer_len}\t{strand_sign}\t.\tID=CRISPR{array_index}_SPACER{s_index}_{s_s}_{s_e};Name=CRISPR{array_index}_REPEAT{s_index}_{s_s}_{s_e};Parent=CRISPR{array_index}_{s_s}_{s_e};Note={r};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                            f.write(line_repeat)
                            f.write(line_spacer)

                        last_repeat = repeats[-1]
                        len_last_repeat = len(last_repeat)
                        last_repeat_start = repeat_starts[-1]
                        last_repeat_end = repeat_ends[-1]
                        last_repeat_index = repeat_indexes[-1]

                        line_repeat = f"{self.acc_num}\tCRISPRidentify\tdirect_repeat\t{last_repeat_start}\t{last_repeat_end}\t{len_last_repeat}\t{strand_sign}\t.\tID=CRISPR{array_index}_REPEAT{last_repeat_index}_{last_repeat_start}_{last_repeat_end};Name=CRISPR{array_index}_REPEAT{array_index}_{last_repeat_start}_{last_repeat_end};Parent=CRISPR{array_index}_{last_repeat_start}_{last_repeat_end};Note={last_repeat};Dbxref=SO:0001459;Ontology_term=CRISPR\n"
                        f.write(line_repeat)


class CompleteFolderSummaryMaker:
    def __init__(self, folder_result):
        self.folder_result = folder_result

        self._list_sub_folders_()
        self._make_complete_summary()

    def _list_sub_folders_(self):
        self.sub_folders = [folder for folder in listdir(self.folder_result)
                            if not isfile(join(self.folder_result, folder))]

    def _make_complete_summary(self):
        summary_path = join(self.folder_result, "Complete_summary.csv")
        global_counter = 1
        flag_header = False
        flag_found_arrays = False
        with open(summary_path, "w") as f:
            for sub_folder_index, sub_folder in enumerate(self.sub_folders):
                complete_path = join(self.folder_result, sub_folder, "Summary.csv")
                if os.path.exists(complete_path):
                    flag_found_arrays = True
                    with open(complete_path) as fr:
                        lines = fr.readlines()
                    if not flag_header:
                        header = lines[0]
                        header = "Name,Global ID," + header
                        f.write(header)
                        flag_header = True

                    for line in lines[1:]:
                        new_line = sub_folder.strip() + "," + str(global_counter) + "," + line
                        global_counter += 1
                        f.write(new_line)

            if not flag_found_arrays:
                f.write("No arrays found")


class CasSummaryMaker:
    def __init__(self, result_path, non_array_data):
        self.result_path = result_path
        self.non_array_data = non_array_data

        self._write_cas_summary()

    def _write_cas_summary(self):
        dict_unstructured_cas = self.non_array_data["Unstructured_Cas"]
        if dict_unstructured_cas:
            result_csv_path = join(self.result_path, 'Cas_Summary.csv')
            with open(result_csv_path, "w") as f:
                f.write(",".join(["ID", "Start", "End", "Length", "Cas type"]))
                f.write("\n")
                counter = 1
                for interval, cas_type in dict_unstructured_cas.items():
                    start, end = interval
                    start, end = int(start), int(end)
                    length = end - start + 1
                    line = ",".join([str(x) for x in [counter, start, end, length, cas_type]]) + "\n"
                    f.write(line)
                    counter += 1


class CompleteCasSummaryFolderMaker:
    def __init__(self, folder_result):
        self.folder_result = folder_result

        self._list_sub_folders_()
        self._make_complete_summary()

    def _list_sub_folders_(self):
        self.sub_folders = [folder for folder in listdir(self.folder_result)
                            if not isfile(join(self.folder_result, folder))]

    def _make_complete_summary(self):
        summary_path = join(self.folder_result, "Complete_Cas_summary.csv")
        global_counter = 1
        flag_header = False
        flag_found_cas = False
        with open(summary_path, "w") as f:
            for sub_folder_index, sub_folder in enumerate(self.sub_folders):
                complete_path = join(self.folder_result, sub_folder, "Cas_Summary.csv")
                if os.path.exists(complete_path):
                    flag_found_cas = True
                    with open(complete_path) as fr:
                        lines = fr.readlines()
                    if not flag_header:
                        header = lines[0]
                        header = "Name,Global ID," + header
                        f.write(header)
                        flag_header = True

                    for line in lines[1:]:
                        new_line = sub_folder.strip() + "," + str(global_counter) + "," + line
                        global_counter += 1
                        f.write(new_line)

            if not flag_found_cas:
                f.write("No cas proteins found")


class FastaOutputArrayMaker:
    def __init__(self, folder_result, categories, non_array_data):
        self.folder_result = folder_result
        self.categories = categories
        self.non_array_data = non_array_data

        self._make_fasta_files()

    def _make_fasta_files(self):
        result_path_repeats = join(self.folder_result, 'Repeats.fasta')
        result_path_spacers = join(self.folder_result, 'Spacers.fasta')
        result_path_arrays = join(self.folder_result, 'Arrays.fasta')
        list_crisprs = [list_info[0][1] for list_info in self.categories[0].values()]

        if list_crisprs:

            with open(result_path_repeats, "w") as fr:
                with open(result_path_spacers, "w") as fs:
                    with open(result_path_arrays, "w") as fa:
                        for index, array in enumerate(list_crisprs):

                            strand = self.non_array_data["Strand"]["Bona-fide"][index]
                            if strand in ("Forward", "Forward (Orientation was not computed)"):
                                output_crispr = OutputCrispr(array)
                            else:
                                output_crispr = OutputCrispr(RevComComputation(array).output())

                            crispr_index = index + 1
                            start, end = output_crispr.start, output_crispr.end

                            repeats = output_crispr.list_repeats
                            spacers = output_crispr.list_spacers
                            array_seq = "".join([r+s for r, s in zip(repeats, spacers)]) + repeats[-1]

                            for repeat_number, repeat in enumerate(repeats, 1):
                                header_line = f">CRISPR_{crispr_index}_{start}_{end}_repeat_{repeat_number}\n"
                                repeat_line = repeat + "\n"
                                fr.write(header_line)
                                fr.write(repeat_line)

                            for spacer_number, spacer in enumerate(spacers, 1):
                                header_line = f">CRISPR_{crispr_index}_{start}_{end}_spacer_{spacer_number}\n"
                                repeat_line = spacer + "\n"
                                fs.write(header_line)
                                fs.write(repeat_line)

                            fa.write(f">CRISPR_{crispr_index}_{start}_{end}\n")
                            fa.write(array_seq + "\n")


class CompleteFastaOutputMaker:
    def __init__(self, folder_result):
        self.folder_result = folder_result

        self._list_sub_folders_()
        self._make_complete_summary()

    def _list_sub_folders_(self):
        self.sub_folders = [folder for folder in listdir(self.folder_result)
                            if not isfile(join(self.folder_result, folder))]

    def _make_complete_summary(self):
        summary_path_repeats = join(self.folder_result, "Complete_repeat_dataset.fasta")
        summary_path_spacers = join(self.folder_result, "Complete_spacer_dataset.fasta")
        summary_path_arrays = join(self.folder_result, "Complete_array_dataset.fasta")

        with open(summary_path_repeats, "w") as fr:
            for sub_folder_index, sub_folder in enumerate(self.sub_folders):
                complete_path = join(self.folder_result, sub_folder, "Repeats.fasta")
                if os.path.exists(complete_path):
                    with open(complete_path) as f:
                        lines = f.readlines()
                        headers = lines[::2]
                        sequences = lines[1::2]
                        for header, sequence in zip(headers, sequences):
                            new_header = ">" + sub_folder.strip() + "_" + header[1:]
                            fr.write(new_header)
                            fr.write(sequence)

        with open(summary_path_spacers, "w") as fs:
            for sub_folder_index, sub_folder in enumerate(self.sub_folders):
                complete_path = join(self.folder_result, sub_folder, "Spacers.fasta")
                if os.path.exists(complete_path):
                    with open(complete_path) as f:
                        lines = f.readlines()
                        headers = lines[::2]
                        sequences = lines[1::2]
                        for header, sequence in zip(headers, sequences):
                            new_header = ">" + sub_folder.strip() + "_" + header[1:]
                            fs.write(new_header)
                            fs.write(sequence)

        with open(summary_path_arrays, "w") as fa:
            for sub_folder_index, sub_folder in enumerate(self.sub_folders):
                complete_path = join(self.folder_result, sub_folder, "Arrays.fasta")
                if os.path.exists(complete_path):
                    with open(complete_path) as f:
                        lines = f.readlines()
                        headers = lines[::2]
                        sequences = lines[1::2]
                        for header, sequence in zip(headers, sequences):
                            new_header = ">" + sub_folder.strip() + "_" + header[1:]
                            fa.write(new_header)
                            fa.write(sequence)










