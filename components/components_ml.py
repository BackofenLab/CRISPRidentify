import operator
import numpy as np
import sklearn
import joblib


class ClassifierWrapper(object):
    def __init__(self, classifier_type, load_option=None, hyper_parameters=None):
        self.classifier_type = classifier_type
        self._hyper_parameters = hyper_parameters
        self._load_option = load_option

        self._init_classifier()

    def _init_classifier(self):
        if self._load_option:
            self._load_model()
        else:
            if self.classifier_type == 'k_near_neighbors':

                if not self._hyper_parameters:
                    self.classifier = sklearn.neighbors.KNeighborsClassifier(n_neighbors=7)
                else:
                    self.classifier = sklearn.neighbors.KNeighborsClassifier(**self._hyper_parameters)

            elif self.classifier_type == 'svm':

                if not self._hyper_parameters:
                    self.classifier = sklearn.svm.SVC()
                else:
                    self.classifier = sklearn.svm.SVC(**self._hyper_parameters)

            elif self.classifier_type == 'naive_bayes':

                if not self._hyper_parameters:
                    self.classifier = sklearn.naive_bayes.GaussianNB()
                else:
                    self.classifier = sklearn.naive_bayes.GaussianNB(**self._hyper_parameters)

            elif self.classifier_type == 'random_forest':

                if not self._hyper_parameters:
                    self.classifier = RandomForestClassifier(max_depth=3, random_state=None)
                else:
                    self.classifier = RandomForestClassifier(**self._hyper_parameters)

            elif self.classifier_type == 'neural_network':

                if not self._hyper_parameters:
                    self.classifier = MLPClassifier(solver='lbfgs', alpha=1e-5,
                                                    hidden_layer_sizes=(100, 100), random_state=None)
                else:
                    self.classifier = MLPClassifier(**self._hyper_parameters)

            elif self.classifier_type == 'extra_trees':

                if not self._hyper_parameters:
                    self.classifier = ExtraTreesClassifier(max_depth=4)
                else:
                    self.classifier = ExtraTreesClassifier(**self._hyper_parameters)

            else:
                raise ValueError('Wrong classifier')

    def _load_model(self):
        self.classifier = joblib.load(self._load_option)

    def train_classifier(self, train_set_pos, train_set_neg):
        train_y_pos = np.ones(len(train_set_pos))
        train_y_neg = np.zeros(len(train_set_neg))
        train_y = np.concatenate([train_y_pos, train_y_neg])
        train_x = np.concatenate([train_set_pos, train_set_neg])
        self.classifier.fit(train_x, train_y)

    def test_classifier(self, test_set_pos, test_set_neg):
        if (test_set_pos is not None) and (test_set_neg is not None):
            test_set_y_pos = np.ones(len(test_set_pos))
            test_set_y_neg = np.zeros(len(test_set_neg))
            test_set_y = np.concatenate([test_set_y_pos, test_set_y_neg])
            test_set_x = np.concatenate([test_set_pos, test_set_neg])

        elif test_set_pos is not None:
            test_set_y = np.ones(len(test_set_pos))
            test_set_x = test_set_pos

        elif test_set_neg is not None:
            test_set_y = np.zeros(len(test_set_neg))
            test_set_x = test_set_neg

        else:
            raise ValueError

        predict = self.classifier.predict(test_set_x)
        dif = test_set_y - predict
        return 1 - np.count_nonzero(dif) / float(len(dif))

    def predict(self, dataset):
        return self.classifier.predict(dataset)

    def predict_proba(self, dataset):
        return self.classifier.predict_proba(dataset)

    def save_model(self, model_name_dot_pkl):
        joblib.dump(self.classifier, model_name_dot_pkl)



