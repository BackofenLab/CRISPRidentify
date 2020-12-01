from components_output_maker import SimpleOutputMaker
from components_output_maker import SummaryOutputMaker
from components_output_maker import SummaryMakerCSV
from components_output_maker import PickleOutputMaker


class OutputMaker:
    def __init__(self, file_path, parameters, result_path, pickle_result_path,
                 categories, non_array_data, list_features, header):
        self.file_path = file_path
        self.parameters = parameters
        self.result_path = result_path
        self.pickle_result_path = pickle_result_path
        self.categories = categories
        self.non_array_data = non_array_data
        self.list_features = list_features
        self.header = header

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

        sm_csv = SummaryMakerCSV(result_path=self.result_path,
                                 categories=self.categories,
                                 non_array_data=self.non_array_data)

        if self.pickle_result_path:
            pom = PickleOutputMaker(file_path=self.file_path,
                                    pickle_result_folder=self.pickle_result_path,
                                    parameters=self.parameters,
                                    categories=self.categories,
                                    non_array_data=self.non_array_data,
                                    header=self.header,
                                    list_feature_names=self.list_features)