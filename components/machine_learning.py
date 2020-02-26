import operator
import numpy as np
import sklearn

from sklearn.externals import joblib
from sklearn import naive_bayes
from sklearn import ensemble
from sklearn import neural_network       
    

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
                
                if not self._hyper_parameters :
                    self.classifier = sklearn.neighbors.KNeighborsClassifier(n_neighbors=7)
                else:
                    self.classifier = sklearn.neighbors.KNeighborsClassifier(**self._hyper_parameters)
                    
            elif self.classifier_type == 'svm':
                
                if not self._hyper_parameters :
                    self.classifier = sklearn.svm.SVC()
                else:
                    self.classifier = sklearn.svm.SVC(**self._hyper_parameters)
                    
            elif self.classifier_type == 'naive_bayes':
                
                if not self._hyper_parameters :
                    self.classifier = sklearn.naive_bayes.GaussianNB()
                else:
                    self.classifier = sklearn.naive_bayes.GaussianNB(**self._hyper_parameters)
                
            elif self.classifier_type == 'random_forest':
                
                if not self._hyper_parameters :
                    self.classifier = RandomForestClassifier(max_depth=3, random_state=None)
                else:
                    self.classifier = RandomForestClassifier(**self._hyper_parameters)
                
            elif self.classifier_type == 'neural_network':
                
                if not self._hyper_parameters :
                    self.classifier = MLPClassifier(solver='lbfgs', alpha=1e-5,
                                             hidden_layer_sizes=(100, 100), random_state=None)
                else:
                    self.classifier = MLPClassifier(**self._hyper_parameters)
                    
            elif self.classifier_type == 'extra_trees':
                
                if not self._hyper_parameters :
                    self.classifier = ExtraTreesClassifier(max_depth=4)
                else:
                    self.classifier = ExtraTreesClassifier(**self._hyper_parameters)
                
            else:
                raise ValueError('Wrong classifier')
    
    def _load_model(self):        
        self.classifier = joblib.load(self._load_option)
        print('You have picked {} as your model'.format(type(self.classifier)))
    
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
        return 1 - np.count_nonzero(dif)/float(len(dif))
    
    def predict(self, dataset):
        return self.classifier.predict(dataset)
    
    def predict_proba(self, dataset):
        return self.classifier.predict_proba(dataset)
    
    def save_model(self, model_name_dot_pkl):
        joblib.dump(self.classifier, model_name_dot_pkl)
        

class AutoML(object):
    def __init__(self, x_train, y_train, x_test, y_test):
        self.train_x = x_train
        self.train_y = y_train
        self.test_x = x_test
        self.test_y = y_test
        
        self.best_models = {}
        
    def _train_k_near_neighbours(self):
        print('training_k_nearest_neighbours')
        best_score = -1
        for neighbours in range(1, 11):
            classifier = sklearn.neighbors.KNeighborsClassifier(n_neighbors=neighbours)
            classifier.fit(self.train_x, self.train_y)
            score = classifier.score(self.test_x, self.test_y)
            if score > best_score:
                best_score = score
                model_info = score, classifier, neighbours
                self.best_models['k_nearest_neighbours'] = model_info        
        
    def _train_svm(self):
        print('training_svm')
        best_score = -1        
        classifier = sklearn.svm.LinearSVC()
        classifier.fit(self.train_x, self.train_y)
        score = classifier.score(self.test_x, self.test_y)
        if score > best_score:
            best_score = score
            model_info = score, classifier, []
            self.best_models['svm'] = model_info
    
    def _train_naive_bayes(self):
        print('training_naive_bayes')
        best_score = -1
        for type_bayes in ['bern', 'gauss', 'multy']:
            if type_bayes == 'bern':
                classifier = sklearn.naive_bayes.BernoulliNB()
            if type_bayes == 'gauss':
                classifier = sklearn.naive_bayes.GaussianNB()
            if type_bayes == 'multy':
                classifier = sklearn.naive_bayes.MultinomialNB()
                
            classifier.fit(self.test_x, self.test_y)
            score = classifier.score(self.train_x, self.train_y)
            if score > best_score:
                best_score = score
                model_info = score, classifier, type_bayes
                self.best_models['naive_bayes'] = model_info            
    
    def _train_random_forest(self):
        print('training_random_forest')
        best_score = -1
        for estimators in [8, 13]:
            for iteration in range(10):
                classifier = sklearn.ensemble.RandomForestClassifier(n_estimators=estimators)
                classifier.fit(self.train_x, self.train_y)
                score = classifier.score(self.test_x, self.test_y)
                if score > best_score:
                    best_score = score
                    model_info = score, classifier, estimators
                    self.best_models['random_forest'] = model_info
                    
    def _train_extra_trees(self):
        print('training_random_forest')
        best_score = -1            
        for iteration in range(10):
            classifier = sklearn.tree.ExtraTreeClassifier()
            classifier.fit(self.test_x, self.test_y)
            score = classifier.score(self.train_x, self.train_y)
            if score > best_score:
                best_score = score
                model_info = score, classifier, []
                self.best_models['extra_trees'] = model_info
    
    def _train_neural_network(self):
        print('training_neural_network')
        best_score = -1
        for hidden_layers in (1, 2):
            if hidden_layers == 1:
                hidden_layer_sizes_options = [(i) for i in range(100, 1100, 100)]
            else:
                hidden_layer_sizes_options = [(i, j) for i in range(100, 1100, 100) for j in range(100, 1100, 100)]
                for hidden_layer_sizes  in hidden_layer_sizes_options:                    
                    classifier = sklearn.neural_network.MLPClassifier(hidden_layer_sizes=hidden_layer_sizes)
                    classifier.fit(self.test_x, self.test_y)
                    score = classifier.score(self.train_x, self.train_y)
                    if score > best_score:
                        best_score = score
                        model_info = score, classifier, hidden_layer_sizes_options
                        self.best_models['neural_network'] = model_info
                        
    
    def train_models(self):
        self._train_k_near_neighbours()
        self._train_svm()
        self._train_naive_bayes()
        self._train_random_forest()
        self._train_extra_trees()
        self._train_neural_network()
        
        print('Done Training')
    
    def pick_best_model(self):
        best_model = max(self.best_models.iteritems(), key=operator.itemgetter(1))[0]
        return best_model, self.best_models[best_model]
        
    

if __name__ == '__main__':
    """some tests"""
    from sklearn.datasets import load_iris
    from sklearn.utils import shuffle

    iris = load_iris()
    X, y = iris.data, iris.target

    print(X.shape)
    print(y.shape)
    indexes = np.where(y<=1)

    X = X[indexes]
    y = y[indexes]

    X, y = shuffle(X, y, random_state=0)

    X_train = X[:90]
    y_train = y[:90]

    X_test = X[90:]
    y_test = y[90:]

    a_ml = AutoML(X_train, y_train, X_test, y_test)
    a_ml.train_models()

    print(a_ml.best_models)
    print('\n'*5)
    print(a_ml.pick_best_model())


