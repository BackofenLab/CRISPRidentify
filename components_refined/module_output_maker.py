from components_output_maker import SimpleOutputMaker
from components_output_maker import SummaryOutputMaker


class OutputMaker:
    def __init__(self, result_path, categories, non_array_data, list_features, header):
        self.result_path = result_path
        self.categories = categories
        self.non_array_data = non_array_data
        self.list_features = list_features
        self.header = header

        self._make_output()

    def _make_output(self):
        som = SimpleOutputMaker(categories=self.categories,
                                result_path=self.result_path,
                                list_features=self.list_features)

        suom = SummaryOutputMaker(result_path=self.result_path,
                                  categories=self.categories,
                                  non_array_data=self.non_array_data,
                                  header=self.header,
                                  list_feature_names=self.list_features)