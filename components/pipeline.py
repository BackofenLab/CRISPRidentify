from module_detection import Detection
from module_detection_refinement import DetectionRefinement
from module_evaluation import ArrayEvaluation
from module_evaluated_arrays_enhancement import EvaluatedArraysEnhancement
from module_non_array_computations import NonArrayComputations
from module_output_maker import OutputMaker


class Pipeline:
    def __init__(self, result_folder_path, pickle_folder_path, file_path, list_ml_classifiers,
                 possible_differentiate_model, flag_possible_differential_model,
                 list_features, parameters, flags):
        self.result_folder_path = result_folder_path + "/" + file_path.split("/")[-1].split(".")[0]
        self.pickle_folder_path = pickle_folder_path
        self.file_path = file_path
        self.list_ml_classifiers = list_ml_classifiers
        self.possible_differentiate_model = possible_differentiate_model
        self.flag_possible_differential_model = flag_possible_differential_model
        self.list_features = [features.strip().split(".") for features in list_features]
        self.flags = flags
        self.parameters = parameters

        self.header = None
        self.dict_fuzzy_crisprs = {}
        self.dict_crispr_candidates = {}
        self.categories = {}
        self.non_array_data = {}

        self._get_header()
        self._run_detection()
        self._run_detection_refinement()
        self._run_evaluation()
        self._results_enhancement()
        self._run_non_crispr_computation()
        self._write_output()

    def _get_header(self):
        with open(self.file_path) as f:
            self.header = f.readline()

    def _run_detection(self):
        print("1. Run initial array detection")
        detection = Detection(file_path=self.file_path,
                              flags=self.flags,
                              parameters=self.parameters)
        self.dict_fuzzy_crisprs = detection.output()

    def _run_detection_refinement(self):
        print("2. Refine detected arrays")
        det_ref = DetectionRefinement(dict_fuzzy_crisprs=self.dict_fuzzy_crisprs,
                                      parameters=self.parameters)
        self.dict_crispr_candidates = det_ref.output()

    def _run_evaluation(self):
        print("3. Evaluate candidates")
        ae = ArrayEvaluation(dict_crispr_array_candidates=self.dict_crispr_candidates,
                             list_ml_classifiers=self.list_ml_classifiers,
                             list_features=self.list_features,
                             parameters=self.parameters)
        self.categories = ae.output()

    def _results_enhancement(self):
        print("4. Enhance evaluated arrays")
        a_enh = EvaluatedArraysEnhancement(file_path=self.file_path,
                                           categories=self.categories,
                                           parameters=self.parameters)
        self.categories = a_enh.output()

    def _run_non_crispr_computation(self):
        print("5. Complement arrays with additional info")
        nac = NonArrayComputations(file_path=self.file_path,
                                   categories=self.categories,
                                   flags_non_arrays_computations=self.flags)
        self.non_array_data = nac.output()

    def _write_output(self):
        print("6. Write down the results")
        om = OutputMaker(file_path=self.file_path,
                         parameters=self.parameters,
                         result_path=self.result_folder_path,
                         pickle_result_path=self.pickle_folder_path,
                         categories=self.categories,
                         non_array_data=self.non_array_data,
                         list_features=self.list_features,
                         header=self.header)






