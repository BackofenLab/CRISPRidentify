from components.components_output_maker import SimpleOutputMaker
from components.components_output_maker import SummaryOutputMaker
from components.components_output_maker import SummaryMakerCSV
from components.components_output_maker import PickleOutputMaker
from components.components_output_maker import CasSummaryMaker
from components.components_output_maker import FastaOutputArrayMaker
from components.components_output_maker import JsonOutputMaker
from components.components_output_maker import GFFOutputMaker
from components.components_output_maker import CompleteFastaOutputMaker
from components.components_output_maker import CompleteFolderSummaryMaker
from components.components_output_maker import CompleteCasSummaryFolderMaker
from components.components_output_maker import SpacerSummaryMaker
from components.components_output_maker import CompleteSpacerCSVMaker


class OutputMaker:
    def __init__(self, file_path, parameters, flags, result_path, pickle_result_path,
                 json_result_path, categories, non_array_data, list_features, header):
        self.file_path = file_path
        self.parameters = parameters
        self.flags = flags
        self.result_path = result_path
        self.pickle_result_path = pickle_result_path
        self.json_result_path = json_result_path
        self.categories = categories
        self.non_array_data = non_array_data
        self.list_features = list_features
        self.header = header
        self.global_result_folder = "/".join(self.result_path.split("/")[:-1])

        self._make_output()

    def _make_output(self):
        som = SimpleOutputMaker(categories=self.categories,
                                result_path=self.result_path,
                                non_array_data=self.non_array_data,
                                list_features=self.list_features)

        suom = SummaryOutputMaker(result_path=self.result_path,
                                  categories=self.categories,
                                  non_array_data=self.non_array_data,
                                  header=self.header,
                                  list_feature_names=self.list_features)

        ssm = SpacerSummaryMaker(categories=self.categories,
                                result_path=self.result_path)

        sm_csv = SummaryMakerCSV(result_path=self.result_path,
                                 categories=self.categories,
                                 non_array_data=self.non_array_data)

        gffom = GFFOutputMaker(result_path=self.result_path, categories=self.categories, non_array_data=self.non_array_data,
                                list_feature_names=self.list_features)
        if self.flags["flag_cas"] is True:
            sm_cas = CasSummaryMaker(result_path=self.result_path,
                                     non_array_data=self.non_array_data)


        #cfsm = CompleteFolderSummaryMaker(folder_result=self.global_result_folder)
        #ccfsm = CompleteCasSummaryFolderMaker(folder_result=self.global_result_folder)

        if self.flags["flag_fasta_report"] is True:
            foam = FastaOutputArrayMaker(folder_result=self.result_path,
                                         categories=self.categories,
                                         non_array_data=self.non_array_data)

            #cfom = CompleteFastaOutputMaker(folder_result=self.global_result_folder)

        if self.pickle_result_path:
            pom = PickleOutputMaker(file_path=self.file_path,
                                    pickle_result_folder=self.pickle_result_path,
                                    parameters=self.parameters,
                                    categories=self.categories,
                                    non_array_data=self.non_array_data,
                                    header=self.header,
                                    list_feature_names=self.list_features)

        if self.json_result_path:
            jom = JsonOutputMaker(file_path=self.file_path,
                                  json_result_folder=self.json_result_path,
                                  categories=self.categories,
                                  non_array_data=self.non_array_data,
                                  list_feature_names=self.non_array_data)
