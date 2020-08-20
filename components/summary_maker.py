from components_rev_com_crispr import RevComComputation


class ResultContainer:
    def __init__(self):
        self.header = None
        self.crispr_results = {}
        self.crispr_rev_com_results = {}

        self.is_element_results = {}

        self.cas_genes_results = {}
        self.dict_cas_genes_separated = {}

        self.strand_results = {}

        self.leader_results = {}
        self.leader_results_rev_com = {}

        self.degenerated_results = {}
        self.degenerated_rev_com_results = {}

        self.merged_results = {}

        self.flag_is_element_calculated = False
        self.flag_cas_genes_calculated = False
        self.flag_strand_calculated = False
        self.flag_degenerated_calculated = False
        self.flag_merged_calculated = False

    def set_crispr_results(self, dict_output_crisprs):
        self.crispr_results = dict_output_crisprs

    def set_crispr_rev_com_results(self, dict_output_crisprs_rev_com):
        self.crispr_rev_com_results = dict_output_crisprs_rev_com

    def set_is_element(self, dict_output_is_elements):
        self.flag_is_element_calculated = True
        self.is_element_results = dict_output_is_elements

    def set_cas_genes(self, dict_cas_genes, dict_cas_genes_separated):
        self.flag_cas_genes_calculated = True
        self.cas_genes_results = dict_cas_genes
        self.dict_cas_genes_separated = dict_cas_genes_separated

    def set_strand(self, dict_strand):
        self.flag_strand_calculated = True
        self.strand_results = dict_strand

    def set_leader(self, dict_leader):
        self.leader_results = dict_leader

    def set_leader_rev_com(self, dict_leader_rev_com):
        self.leader_results_rev_com = dict_leader_rev_com

    def set_degenerated(self, dict_degenerated):
        self.flag_degenerated_calculated = True
        self.degenerated_results = dict_degenerated

    def set_degenerated_rev_com(self, dict_degenerated_rev_com):
        self.degenerated_rev_com_results = dict_degenerated_rev_com

    def set_merged(self, dict_merged):
        self.flag_merged_calculated = True
        self.merged_results = {key: OutputCrispr(value) for key, value in dict_merged.items()}


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


class SummaryMaker:
    def __init__(self, header, list_crisprs, list_scores, list_feature_names, list_feature_vectors,
                 additional_calculations, flag_is, flag_cas, flag_degenerated, acc_num):
        self.acc_num = acc_num
        self.header = header
        self.list_scores = list_scores
        self.list_feature_vectors = list_feature_vectors
        self.list_feature_names = list_feature_names
        self.list_crisprs = list_crisprs
        self.additional_calculations = additional_calculations
        self.flag_is = flag_is
        self.flag_cas = flag_cas
        self.flag_degenerated = flag_degenerated

        self.result_container = ResultContainer()

        self._perform_calculations()

    def _perform_calculations(self):
        dict_crisprs = {index: OutputCrispr(crispr_candidate=crispr_candidate)
                        for index, crispr_candidate in enumerate(self.list_crisprs)}

        dict_crisprs_rev_com = {index: OutputCrispr(crispr_candidate=RevComComputation(crispr_candidate).output())
                                for index, crispr_candidate in enumerate(self.list_crisprs)}

        self.result_container.set_crispr_results(dict_crisprs)
        self.result_container.set_crispr_rev_com_results(dict_crisprs_rev_com)

        if self.flag_is:
            is_data = self.additional_calculations.calculate_is_elements()
            self.result_container.set_is_element(is_data)

        if self.flag_cas:
            cas_data, cas_data_separated = self.additional_calculations.calculate_cas_proteins()
            self.result_container.set_cas_genes(cas_data, cas_data_separated)

        if self.flag_degenerated:
            degenerated_original = self.additional_calculations.calculate_degenerated()
            degenerated_as_output = {key: OutputCrispr(value) for key, value in degenerated_original.items()}
            degenerated_rev_com = {key: OutputCrispr(RevComComputation(value).output())
                                   for key, value in degenerated_original.items()}

            self.result_container.set_degenerated(degenerated_as_output)
            self.result_container.set_degenerated_rev_com(degenerated_rev_com)

        strand_data = self.additional_calculations.calculate_strand()
        self.result_container.set_strand(strand_data)

        leader_data = self.additional_calculations.calculate_leader()
        self.result_container.set_leader(leader_data)

        leader_data_rev_com = self.additional_calculations.calculate_leader_rev_com()
        self.result_container.set_leader_rev_com(leader_data_rev_com)

        merged_crisprs = self.additional_calculations.calculate_merged()
        self.result_container.set_merged(merged_crisprs)

    def write_text_summary(self, file_name):
        with open(file_name, "w") as f:
            f.write(self.header)
            f.write("\n")

            if self.flag_degenerated:
                crispr_container = self.result_container.degenerated_results
                crispr_container_rev_com = self.result_container.degenerated_rev_com_results
            else:
                crispr_container = self.result_container.crispr_results
                crispr_container_rev_com = self.result_container.crispr_rev_com_results

            for index in sorted(crispr_container.keys()):

                if index in self.result_container.cas_genes_results:
                    cas_genes = self.result_container.cas_genes_results[index]
                    if cas_genes:
                        f.write("Cas genes: ")
                        for index_cas, cluster in enumerate(cas_genes):
                            if index_cas != (len(cas_genes) - 1):
                                f.write("{} [{}-{}], ".format(cluster[2], cluster[0], cluster[1]))
                            else:
                                f.write("{} [{}-{}]".format(cluster[2], cluster[0], cluster[1]))
                    f.write("\n\n")

                strand = self.result_container.strand_results[index]
                if strand == "Reverse":
                    output_crispr = crispr_container_rev_com[index]
                else:
                    output_crispr = crispr_container[index]

                crispr_index = index + 1
                start, end = output_crispr.start, output_crispr.end
                number_of_repeats = output_crispr.number_of_repeats
                avg_length_repeat, avg_length_spacer = output_crispr.avg_repeat_length, output_crispr.avg_spacer_length

                if index > 0:
                    f.write('\n{}\n\n'.format('=' * 100))
                f.write("CRISPR: {}, {}-{}, number of Repeats: {}, avg. length of Repeat: {}, avg length of Spacer: {}\n\n"
                        .format(str(crispr_index), start, end, number_of_repeats,
                                avg_length_repeat, avg_length_spacer))

                f.write(output_crispr.dot_representation)

                f.write("\n")

                f.write("Leader region\n")

                if strand == "Reverse":
                    leader = self.result_container.leader_results_rev_com[index]
                else:
                    leader = self.result_container.leader_results[index]

                f.write(leader)

                f.write("\n\nStrand: {}\n\n".format(self.result_container.strand_results[index]))

                f.write("#   Array features:\n")
                list_reported_features = []
                for group_of_features, feature_vector in zip(self.list_feature_names, self.list_feature_vectors[index]):
                    for feature_name, feature_value in zip(group_of_features, feature_vector[0]):
                        if feature_name not in list_reported_features:
                            f.write("#   {}: {}\n".format(feature_name, feature_value))
                            list_reported_features.append(feature_name)

                f.write("_" * 30)
                f.write("\n")
                f.write("#   Certainty Score: {}\n\n".format(self.list_scores[index]))

                if index in self.result_container.is_element_results:
                    f.write("IS Element: {} [{}-{}]\n\n".format(self.result_container.is_element_results[index][4],
                                                                 self.result_container.is_element_results[index][0],
                                                                 self.result_container.is_element_results[index][1]))

            last_index = len(self.list_crisprs)
            if last_index in self.result_container.cas_genes_results:
                f.write("\n\n")
                cas_genes = self.result_container.cas_genes_results[last_index]
                if cas_genes:
                    f.write("Cas genes: ")
                    for index, cluster in enumerate(cas_genes):
                        if index != (len(cas_genes) - 1):
                            f.write("{} [{}-{}], ".format(cluster[2], cluster[0], cluster[1]))
                        else:
                            f.write("{} [{}-{}]".format(cluster[2], cluster[0], cluster[1]))
                f.write("\n\n")

    def write_json_summary(self, file_name):
        pass

    def write_merged_summary(self, file_name):
        with open(file_name, "w") as f:
            for index, cc in self.result_container.merged_results.items():
                f.write(cc.dot_representation)
                f.write("\n\n\n")

    def write_markdown_summary(self, filename):
        pass

    def write_bed_summary(self, filename):
        with open(filename, "w") as f:
            if self.flag_degenerated:
                crispr_container = self.result_container.degenerated_results
            else:
                crispr_container = self.result_container.crispr_results
            for index in sorted(crispr_container.keys()):
                output_crispr = crispr_container[index]

                if index in self.result_container.dict_cas_genes_separated:
                    clusters = self.result_container.dict_cas_genes_separated[index]
                    for cluster in clusters:
                        for cas in cluster:
                            f.write("{}\t{}\t{}\t{}\t.\t{}\n".format(self.acc_num, cas[0], cas[1], cas[2], cas[3]))

                start, end = output_crispr.start, output_crispr.end
                strand = self.result_container.strand_results[index]
                f.write("{}\t{}\t{}\t{}\t.\t{}\n".format(self.acc_num, start, end, "CRISPR" + str(index + 1), strand))

                if index in self.result_container.is_element_results:
                    is_start = self.result_container.is_element_results[index][0]
                    is_end = self.result_container.is_element_results[index][1]
                    is_strand = self.result_container.is_element_results[index][1]
                    f.write("{}\t{}\t{}\t{}\t.\t{}\n".format(self.acc_num, is_start, is_end, "IS", is_strand))

            last_index = len(self.list_crisprs)
            if last_index in self.result_container.dict_cas_genes_separated:
                clusters = self.result_container.dict_cas_genes_separated[last_index]
                for cluster in clusters:
                    for cas in cluster:
                        f.write("{}\t{}\t{}\t{}\t.\t{}\n".format(self.acc_num, cas[0], cas[1], cas[2], cas[3]))

    def write_gff_summary(self, filename):
        with open(filename, "w") as f:
            if self.flag_degenerated:
                crispr_container = self.result_container.degenerated_results
            else:
                crispr_container = self.result_container.crispr_results
            for index in sorted(crispr_container.keys()):
                output_crispr = crispr_container[index]
                main_line_ending = [f"ID=CRISPR{index + 1}_{output_crispr.start}_{output_crispr.end}",
                                    f"Note={output_crispr.consensus}",
                                    "DBxref=XXXXXX",
                                    "Ontology_term=CRISPR"]
                main_line_ending = ";".join(main_line_ending)
                main_line = [self.acc_num,
                             "CRISPRIdentifier",
                             "repeat_region",
                             output_crispr.start,
                             output_crispr.end,
                             output_crispr.end - output_crispr.start,
                             "+", ".",
                             main_line_ending]
                main_line = "\t".join(main_line)
                f.write(main_line)
                f.write("\n")

                for repeat_index in range(len(output_crispr.list_repeats)):
                    repeat = output_crispr.list_repeats[repeat_index]
                    spacer = output_crispr.list_spacers[repeat_index]

                    repeat_start = output_crispr.list_repeat_starts[repeat_index]
                    repeat_end = repeat_start + len(repeat)

                    spacer_start = output_crispr.list_spacer_starts[repeat_index]
                    spacer_end = spacer_start + len(spacer)

                    repeat_line_ending = [f"ID=CRISPR{index+1}_REPEAT{repeat_index+1}_{repeat_start}_{repeat_end}",
                                          f"Name=CRISPR{index+1}_REPEAT{repeat_index+1}_{repeat_start}_{repeat_end}",
                                          f"Parent=CRISPR{index + 1}_{output_crispr.start}_{output_crispr.end}"]

                    repeat_line_ending = ";".join(repeat_line_ending)
                    repeat_line = [self.acc_num,
                                   "CRISPRIdentifier",
                                   "direct_repeat",
                                   repeat_start,
                                   repeat_end,
                                   repeat_end - repeat_start,
                                   "+", ".",
                                   repeat_line_ending]

                    repeat_line = "\t".join(repeat_line)

                    f.write(repeat_line)
                    f.write("\n")

                    spacer_line_ending = [f"ID=CRISPR{index + 1}_SPACER{repeat_index + 1}_{repeat_start}_{repeat_end}",
                                          f"Name=CRISPR{index + 1}_SPACER{repeat_index + 1}_{repeat_start}_{repeat_end}",
                                          f"Parent=CRISPR{index + 1}_{output_crispr.start}_{output_crispr.end}"]

                    spacer_line = [self.acc_num,
                                   "CRISPRIdentifier",
                                   "spacer",
                                   spacer_start,
                                   spacer_end,
                                   spacer_end - spacer_start,
                                   "+", ".",
                                   spacer_line_ending]

                    f.write(spacer_line)
                    f.write("\n")

















