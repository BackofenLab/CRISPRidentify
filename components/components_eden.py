from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import collections


__magic__ = ('EDeN', 42)
__magic_py2hash__ = -7048895691955021301
__magic_py3hash__ = -1821860980875793120

_bitmask_ = 4294967295

import joblib
import networkx as nx
import math
import numpy as np
from sklearn import metrics
from sklearn.cluster import MiniBatchKMeans
from scipy import stats
from scipy.sparse import csr_matrix
from scipy.sparse import vstack
from collections import defaultdict, deque
import dill
from sklearn.base import BaseEstimator, TransformerMixin
from scipy import io
import requests
import os
import sys
from collections import deque
from itertools import tee
import random
import logging.handlers

import multiprocessing as mp
import time

from itertools import tee
import logging
logger = logging.getLogger(__name__)


class AbstractVectorizer(BaseEstimator, TransformerMixin):
    """Interface declaration for the Vectorizer class."""

    def annotate(self, graphs, estimator=None, reweight=1.0, relabel=False):
        raise NotImplementedError("Should have implemented this")

    def set_params(self, **args):
        raise NotImplementedError("Should have implemented this")

    def transform(self, graphs):
        raise NotImplementedError("Should have implemented this")

    def vertex_transform(self, graph):
        raise NotImplementedError("Should have implemented this")


def run_dill_encoded(what):
    """Use dill as replacement for pickle to enable multiprocessing on instance methods"""

    fun, args = dill.loads(what)
    return fun(*args)


def apply_async(pool, fun, args, callback=None):
    """
    Wrapper around apply_async() from multiprocessing, to use dill instead of pickle.
    This is a workaround to enable multiprocessing of classes.
    """
    return pool.apply_async(run_dill_encoded, (dill.dumps((fun, args)),), callback=callback)


def fast_hash_2(dat_1, dat_2, bitmask=_bitmask_):
    return int(hash((dat_1, dat_2)) & bitmask) + 1


def fast_hash_3(dat_1, dat_2, dat_3, bitmask=_bitmask_):
    return int(hash((dat_1, dat_2, dat_3)) & bitmask) + 1


def fast_hash_4(dat_1, dat_2, dat_3, dat_4, bitmask=_bitmask_):
    return int(hash((dat_1, dat_2, dat_3, dat_4)) & bitmask) + 1


def fast_hash(vec, bitmask=_bitmask_):
    return int(hash(tuple(vec)) & bitmask) + 1


def fast_hash_vec(vec, bitmask=_bitmask_):
    hash_vec = []
    running_hash = 0xAAAAAAAA
    for i, vec_item in enumerate(vec):
        running_hash ^= hash((running_hash, vec_item, i))
        hash_vec.append(int(running_hash & bitmask) + 1)
    return hash_vec


def auto_label(graphs, n_clusters=16, **opts):
    """Label nodes with cluster id.

    Cluster nodes using as features the output of vertex_vectorize.
    """
    data_list = Vectorizer(**opts).vertex_transform(graphs)
    data_matrix = vstack(data_list)
    clu = MiniBatchKMeans(n_clusters=n_clusters, n_init=10)
    clu.fit(data_matrix)
    preds = clu.predict(data_matrix)
    vecs = clu.transform(data_matrix)
    sizes = [m.shape[0] for m in data_list]
    label_list = []
    vecs_list = []
    pointer = 0
    for size in sizes:
        label_list.append(preds[pointer: pointer + size])
        vecs_list.append(vecs[pointer: pointer + size])
        pointer += size
    return label_list, vecs_list


def auto_relabel(graphs, n_clusters=16, **opts):
    """Label nodes with cluster id."""
    graphs, graphs_ = tee(graphs)
    label_list, vecs_list = auto_label(graphs_, n_clusters=n_clusters, **opts)
    relabeled_graphs = []
    for labels, vecs, orig_graph in zip(label_list, vecs_list, graphs):
        graph = nx.Graph(orig_graph)
        for label, vec, u in zip(labels, vecs, graph.nodes()):
            graph.nodes[u]['label'] = label
            graph.nodes[u]['vec'] = list(vec)
        relabeled_graphs.append(graph)
    return relabeled_graphs


def vectorize(graphs, **opts):
    """Transform real vector labeled, weighted graphs in sparse vectors."""
    return Vectorizer(**opts).transform(graphs)


def vertex_vectorize(graphs, **opts):
    """Transform a list of networkx graphs into a list of sparse matrices."""
    return Vectorizer(**opts).vertex_transform(graphs)


def annotate(graphs,
             estimator=None,
             reweight=1.0,
             vertex_features=False, **opts):
    """Return graphs with extra node attributes: importance and features."""
    result = Vectorizer(**opts).annotate(graphs,
                                         estimator=estimator,
                                         reweight=reweight,
                                         vertex_features=vertex_features)
    return list(result)


def kernel_matrix(graphs, **opts):
    """Return the kernel matrix."""
    data_matrix = vectorize(graphs, **opts)
    return metrics.pairwise.pairwise_kernels(data_matrix, metric='linear')

# --------------------------------------------------------------------------


class Vectorizer(AbstractVectorizer):
    """Transform real vector labeled, weighted graphs in sparse vectors."""

    def __init__(self,
                 complexity=3,
                 r=None,
                 d=None,
                 min_r=0,
                 min_d=0,
                 weights_dict=None,
                 auto_weights=False,
                 nbits=16,
                 normalization=True,
                 inner_normalization=True,
                 positional=False,
                 discrete=True,
                 key_label='label',
                 key_weight='weight',
                 key_nesting='nesting',
                 key_importance='importance',
                 key_class='class',
                 key_vec='vec',
                 key_svec='svec'):
        """Constructor.

        Parameters
        ----------
        complexity : int (default 3)
            The complexity of the features extracted.
            This is equivalent to setting r = complexity, d = complexity.

        r : int
            The maximal radius size.

        d : int
            The maximal distance size.

        min_r : int
            The minimal radius size.

        min_d : int
            The minimal distance size.

        weights_dict : dict of floats
            Dictionary with keys = pairs (radius, distance) and
            value = weights.

        auto_weights : bool (default False)
            Flag to set to 1 the weight of the kernels for r=i, d=i
            for i in range(complexity)

        nbits : int (default 16)
            The number of bits that defines the feature space size:
            |feature space|=2^nbits.

        normalization : bool (default True)
            Flag to set the resulting feature vector to have unit euclidean
            norm.

        inner_normalization : bool (default True)
            Flag to set the feature vector for a specific combination of the
            radius and distance size to have unit euclidean norm.
            When used together with the 'normalization' flag it will be applied
            first and then the resulting feature vector will be normalized.

        positional : bool (default False)
            Flag to make the relative position be sorted by the node ID value.
            This is useful for ensuring isomorphism for sequences.

        discrete: bool (default False)
            Flag to activate more efficient computation of vectorization
            considering only discrete labels and ignoring vector attributes.

        key_label : string (default 'label')
            The key used to indicate the label information in nodes.

        key_weight : string (default 'weight')
            The key used to indicate the weight information in nodes.

        key_nesting : string (default 'nesting')
            The key used to indicate the nesting type in edges.

        key_importance : string (default 'importance')
            The key used to indicate the importance information in nodes.

        key_class : string (default 'class')
            The key used to indicate the predicted class associated to
            the node.

        key_vec : string (default 'vec')
            The key used to indicate the vector label information in nodes.

        key_svec : string (default 'svec')
            The key used to indicate the sparse vector label information
            in nodes.
        """
        self.name = self.__class__.__name__
        self.__version__ = '1.0.1'
        self.complexity = complexity
        if r is None:
            r = complexity
        if d is None:
            d = complexity
        self.r = r
        self.d = d
        self.min_r = min_r
        self.min_d = min_d
        self.weights_dict = weights_dict
        if auto_weights:
            self.weights_dict = {(i, i + 1): 1 for i in range(max(r, d))}
        self.nbits = nbits
        self.normalization = normalization
        self.inner_normalization = inner_normalization
        self.positional = positional
        self.discrete = discrete
        self.bitmask = pow(2, nbits) - 1
        self.feature_size = self.bitmask + 2
        self.key_label = key_label
        self.key_weight = key_weight
        self.key_nesting = key_nesting
        self.key_importance = key_importance
        self.key_class = key_class
        self.key_vec = key_vec
        self.key_svec = key_svec

    def set_params(self, **args):
        """Set the parameters of the vectorizer."""
        if args.get('complexity', None) is not None:
            self.complexity = args['complexity']
            self.r = self.complexity
            self.d = self.complexity
        if args.get('r', None) is not None:
            self.r = args['r']
        if args.get('d', None) is not None:
            self.d = args['d']
        if args.get('min_r', None) is not None:
            self.min_r = args['min_r']
        if args.get('min_d', None) is not None:
            self.min_d = args['min_d']
        if args.get('nbits', None) is not None:
            self.nbits = args['nbits']
            self.bitmask = pow(2, self.nbits) - 1
            self.feature_size = self.bitmask + 2
        if args.get('normalization', None) is not None:
            self.normalization = args['normalization']
        if args.get('inner_normalization', None) is not None:
            self.inner_normalization = args['inner_normalization']
        if args.get('positional', None) is not None:
            self.positional = args['positional']

    def get_params(self):
        """Get parameters for teh vectorizer.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        return self.__dict__

    def __repr__(self):
        """string."""
        return serialize_dict(self.__dict__, offset='large')

    def save(self, model_name):
        """save."""
        joblib.dump(self, model_name, compress=1)

    def load(self, obj):
        """load."""
        self.__dict__.update(joblib.load(obj).__dict__)

    def transform(self, graphs):
        """Transform a list of networkx graphs into a sparse matrix.

        Parameters
        ----------
        graphs : list[graphs]
            The input list of networkx graphs.

        Returns
        -------
        data_matrix : array-like, shape = [n_samples, n_features]
            Vector representation of input graphs.

        >>> # transforming the same graph
        >>> import networkx as nx
        >>> def get_path_graph(length=4):
        ...     g = nx.path_graph(length)
        ...     for n,d in g.nodes(data=True):
        ...         d['label'] = 'C'
        ...     for a,b,d in g.edges(data=True):
        ...         d['label'] = '1'
        ...     return g
        >>> g = get_path_graph(4)
        >>> g2 = get_path_graph(5)
        >>> g2.remove_node(4)
        >>> v = Vectorizer()
        >>> def vec_to_hash(vec):
        ...     return hash(tuple(vec.data + vec.indices))
        >>> vec_to_hash(v.transform([g])) == vec_to_hash(v.transform([g2]))
        True
        """
        instance_id = None
        feature_rows = []
        for instance_id, graph in enumerate(graphs):
            self._test_goodness(graph)
            feature_rows.append(self._transform(graph))
        if instance_id is None:
            raise Exception('ERROR: something went wrong:\
                no graphs are present in current iterator.')
        data_matrix = self._convert_dict_to_sparse_matrix(feature_rows)
        return data_matrix

    def vertex_transform(self, graphs):
        """Transform a list of networkx graphs into a list of sparse matrices.

        Each matrix has dimension n_nodes x n_features, i.e. each vertex is
        associated to a sparse vector that encodes the neighborhood of the
        vertex up to radius + distance.

        Parameters
        ----------
        graphs : list[graphs]
            The input list of networkx graphs.

        Returns
        -------
        matrix_list : array-like, shape = [n_samples, [n_nodes, n_features]]
            Vector representation of each vertex in the input graphs.

        """
        matrix_list = []
        for instance_id, graph in enumerate(graphs):
            self._test_goodness(graph)
            graph = self._graph_preprocessing(graph)
            # extract per vertex feature representation
            data_matrix = self._compute_vertex_based_features(graph)
            matrix_list.append(data_matrix)
        return matrix_list

    def _test_goodness(self, graph):
        if graph.number_of_nodes() == 0:
            raise Exception('ERROR: something went wrong, empty graph.')

    def _convert_dict_to_sparse_matrix(self, feature_rows):
        if len(feature_rows) == 0:
            raise Exception('ERROR: something went wrong, empty features.')
        data, row, col = [], [], []
        for i, feature_row in enumerate(feature_rows):
            if len(feature_row) == 0:
                # case of empty feature set for a specific instance
                row.append(i)
                col.append(0)
                data.append(0)
            else:
                for feature in feature_row:
                    row.append(i)
                    col.append(feature)
                    data.append(feature_row[feature])
        shape = (max(row) + 1, self.feature_size)
        data_matrix = csr_matrix((data, (row, col)),
                                 shape=shape, dtype=np.float64)
        return data_matrix

    def _init_weight_preprocessing(self, graph):
        graph.graph['weighted'] = False
        for u in graph.nodes():
            if graph.nodes[u].get(self.key_weight, False):
                graph.graph['weighted'] = True
                break

    def _weight_preprocessing(self, graph):
        # if at least one vertex or edge is weighted then ensure that all
        # vertices and edges are weighted in this case use a default weight
        # of 1 if the weight attribute is missing
        self._init_weight_preprocessing(graph)
        if graph.graph['weighted'] is True:
            for u in graph.nodes():
                if graph.nodes[u].get(self.key_weight, False) is False:
                    graph.nodes[u][self.key_weight] = 1

    def _graph_preprocessing(self, original_graph):
        graph = _edge_to_vertex_transform(original_graph)
        self._weight_preprocessing(graph)
        _label_preprocessing(graph,
                             key_label=self.key_label,
                             bitmask=self.bitmask)
        self._compute_distant_neighbours(graph, max(self.r, self.d) * 2)
        self._compute_neighborhood_graph_hash_cache(graph)
        if graph.graph.get('weighted', False):
            self._compute_neighborhood_graph_weight_cache(graph)
        return graph

    def _transform(self, original_graph):
        graph = self._graph_preprocessing(original_graph)
        # collect all features for all vertices for each label_index
        feature_list = defaultdict(lambda: defaultdict(float))
        for v in graph.nodes():
            # only for vertices of type 'node', i.e. not for the 'edge' type
            if graph.nodes[v].get('node', False):
                self._transform_vertex(graph, v, feature_list)
        _clean_graph(graph)
        return self._normalization(feature_list)

    def _update_feature_list(self, node_feature_list, feature_list):
        for radius_dist_key in node_feature_list:
            for feature in node_feature_list[radius_dist_key]:
                val = node_feature_list[radius_dist_key][feature]
                feature_list[radius_dist_key][feature] += val

    def _transform_vertex(self, graph, vertex_v, feature_list):
        if self.discrete:
            # for all distances
            root_dist_dict = graph.nodes[vertex_v]['remote_neighbours']
            for distance in range(self.min_d * 2, (self.d + 1) * 2, 2):
                if distance in root_dist_dict:
                    node_set = root_dist_dict[distance]
                    for vertex_u in node_set:
                        if graph.nodes[vertex_u].get('node', False):
                            self._transform_vertex_pair(
                                graph, vertex_v, vertex_u,
                                distance, feature_list)
            self._transform_vertex_nesting(graph, vertex_v, feature_list)
        else:
            node_feature_list = defaultdict(lambda: defaultdict(float))
            # for all distances
            root_dist_dict = graph.nodes[vertex_v]['remote_neighbours']
            for distance in range(self.min_d * 2, (self.d + 1) * 2, 2):
                if distance in root_dist_dict:
                    node_set = root_dist_dict[distance]
                    for vertex_u in node_set:
                        if graph.nodes[vertex_u].get('node', False):
                            self._transform_vertex_pair(
                                graph,
                                vertex_v,
                                vertex_u,
                                distance,
                                node_feature_list)
            self._transform_vertex_nesting(
                graph,
                vertex_v,
                node_feature_list)
            node_feature_list = self._add_vector_labes(
                graph,
                vertex_v,
                node_feature_list)
            self._update_feature_list(
                node_feature_list,
                feature_list)
            node_sparse_feature_list = self._add_sparse_vector_labes(
                graph,
                vertex_v,
                node_feature_list)
            self._update_feature_list(
                node_sparse_feature_list,
                feature_list)

    def _add_vector_labes(self, graph, vertex_v, node_feature_list):
        # add the vector with an offset given by the feature, multiplied by val
        vec = graph.nodes[vertex_v].get(self.key_vec, None)
        if vec:
            vec_feature_list = defaultdict(lambda: defaultdict(float))
            for radius_dist_key in node_feature_list:
                for feature in node_feature_list[radius_dist_key]:
                    val = node_feature_list[radius_dist_key][feature]
                    for i, vec_val in enumerate(vec):
                        key = (feature + i) % self.bitmask
                        vec_feature_list[radius_dist_key][key] += val * vec_val
            node_feature_list = vec_feature_list
        return node_feature_list

    def _add_sparse_vector_labes(self, graph, vertex_v, node_feature_list):
        # add the vector with a feature resulting from hashing
        # the discrete labeled graph sparse encoding with the sparse vector
        # feature, the val is then multiplied.
        svec = graph.nodes[vertex_v].get(self.key_svec, None)
        if svec:
            vec_feature_list = defaultdict(lambda: defaultdict(float))
            for radius_dist_key in node_feature_list:
                for feature in node_feature_list[radius_dist_key]:
                    val = node_feature_list[radius_dist_key][feature]
                    for i in svec:
                        vec_val = svec[i]
                        key = fast_hash_2(feature, i, self.bitmask)
                        vec_feature_list[radius_dist_key][key] += val * vec_val
            node_feature_list = vec_feature_list
        return node_feature_list

    def _transform_vertex_nesting(self, graph, vertex_v, feature_list):
        # find all vertices, if any, that are second point of nesting edge
        endpoints = self._find_second_endpoint_of_nesting_edge(graph, vertex_v)
        for endpoint, connection_weight in endpoints:
            # for all vertices at distance d from each such second endpoint
            endpoint_dist_dict = graph.nodes[endpoint]['remote_neighbours']
            for distance in range(self.min_d * 2, (self.d + 1) * 2, 2):
                if distance in endpoint_dist_dict:
                    node_set = endpoint_dist_dict[distance]
                    # for all nodes u at distance distance from endpoint
                    for vertex_u in node_set:
                        if graph.nodes[vertex_u].get('node', False):
                            self._transform_vertex_pair(
                                graph, vertex_v, vertex_u,
                                distance, feature_list,
                                connection_weight=connection_weight)

    def _find_second_endpoint_of_nesting_edge(self, graph, vertex_v):
        endpoints = []
        # find all neighbors
        for u in graph.neighbors(vertex_v):
            # test for type
            if graph.nodes[u].get(self.key_nesting, False):
                # if type is nesting
                # find endpoint that is not original vertex_v
                vertices = [j for j in graph.neighbors(u) if j != vertex_v]
                assert(len(vertices) == 1)
                connection_weight = graph.nodes[u].get(self.key_weight, 1)
                endpoints.append((vertices[0], connection_weight))
        return endpoints

    def _transform_vertex_pair(self,
                               graph,
                               vertex_v,
                               vertex_u,
                               distance,
                               feature_list,
                               connection_weight=1):
        cw = connection_weight
        # for all radii
        for radius in range(self.min_r * 2, (self.r + 1) * 2, 2):
            # Note: to be compatible with external radius, distance
            # we need to revert to r/2 and d/2
            radius_dist_key = (radius / 2, distance / 2)
            if self.weights_dict is None or \
                    self.weights_dict.get(radius_dist_key, 0) != 0:
                self._transform_vertex_pair_valid(graph,
                                                  vertex_v,
                                                  vertex_u,
                                                  radius,
                                                  distance,
                                                  feature_list,
                                                  connection_weight=cw)

    def _transform_vertex_pair_valid(self,
                                     graph,
                                     vertex_v,
                                     vertex_u,
                                     radius,
                                     distance,
                                     feature_list,
                                     connection_weight=1):
        cw = connection_weight
        # we need to revert to r/2 and d/2
        radius_dist_key = (radius / 2, distance / 2)
        # reweight using external weight dictionary
        len_v = len(graph.nodes[vertex_v]['neigh_graph_hash'])
        len_u = len(graph.nodes[vertex_u]['neigh_graph_hash'])
        if radius < len_v and radius < len_u:
            # feature as a pair of neighborhoods at a radius,distance
            # canonicalization of pair of neighborhoods
            vertex_v_labels = graph.nodes[vertex_v]['neigh_graph_hash']
            vertex_v_hash = vertex_v_labels[radius]
            vertex_u_labels = graph.nodes[vertex_u]['neigh_graph_hash']
            vertex_u_hash = vertex_u_labels[radius]
            if vertex_v_hash < vertex_u_hash:
                first_hash, second_hash = (vertex_v_hash, vertex_u_hash)
            else:
                first_hash, second_hash = (vertex_u_hash, vertex_v_hash)
            feature = fast_hash_4(
                first_hash, second_hash, radius, distance, self.bitmask)
            # half features are those that ignore the central vertex v
            # the reason to have those is to help model the context
            # independently from the identity of the vertex itself
            half_feature = fast_hash_3(vertex_u_hash,
                                       radius, distance, self.bitmask)
            if graph.graph.get('weighted', False) is False:
                feature_list[radius_dist_key][feature] += cw
                feature_list[radius_dist_key][half_feature] += cw
            else:
                weight_v = graph.nodes[vertex_v]['neigh_graph_weight']
                weight_u = graph.nodes[vertex_u]['neigh_graph_weight']
                weight_vu_radius = weight_v[radius] + weight_u[radius]
                val = cw * weight_vu_radius
                # Note: add a feature only if the value is not 0
                if val != 0:
                    feature_list[radius_dist_key][feature] += val
                    half_val = cw * weight_u[radius]
                    feature_list[radius_dist_key][half_feature] += half_val

    def _normalization(self, feature_list):
        # inner normalization per radius-distance
        feature_vector = {}
        for r_d_key in feature_list:
            features = feature_list[r_d_key]
            norm = 0
            for count in features.values():
                norm += count * count
            sqrt_norm = math.sqrt(norm)
            if self.weights_dict is not None:
                # reweight using external weight dictionary
                if self.weights_dict.get(r_d_key, None) is not None:
                    sqrtw = math.sqrt(self.weights_dict[r_d_key])
                    sqrt_norm = sqrt_norm / sqrtw
            for feature_id, count in features.items():
                if self.inner_normalization:
                    feature_vector_value = float(count) / sqrt_norm
                else:
                    feature_vector_value = count
                feature_vector[feature_id] = feature_vector_value
        # global normalization
        if self.normalization:
            normalized_feature_vector = {}
            total_norm = 0.0
            for value in feature_vector.values():
                total_norm += value * value
            sqrt_total_norm = math.sqrt(float(total_norm))
            for feature_id, value in list(feature_vector.items()):
                feature_vector_value = value / sqrt_total_norm
                normalized_feature_vector[feature_id] = feature_vector_value
            return normalized_feature_vector
        else:
            return feature_vector

    def _compute_neighborhood_graph_hash_cache(self, graph):
        assert (len(graph) > 0), 'ERROR: Empty graph'
        for u in graph.nodes():
            if graph.nodes[u].get('node', False):
                self._compute_neighborhood_graph_hash(u, graph)

    def _compute_neighborhood_graph_hash(self, root, graph):
        # list all hashed labels at increasing distances
        hash_list = []
        # for all distances
        root_dist_dict = graph.nodes[root]['remote_neighbours']
        for node_set in root_dist_dict.keys():
            # create a list of hashed labels
            hash_label_list = []
            for v in root_dist_dict[node_set]:
                # compute the vertex hashed label by hashing the hlabel
                # field
                # with the degree of the vertex (obtained as the size of
                # the adjacency dictionary for the vertex v)
                # or, in case positional is set, using the relative
                # position of the vertex v w.r.t. the root vertex
                if self.positional:
                    vhlabel = fast_hash_2(
                        graph.nodes[v]['hlabel'], root - v)
                else:
                    vhlabel = fast_hash_2(
                        graph.nodes[v]['hlabel'], len(graph[v]))
                hash_label_list.append(vhlabel)
            # sort it
            hash_label_list.sort()
            # hash it
            hashed_nodes_at_distance_d_in_neighborhood = fast_hash(
                hash_label_list)
            hash_list.append(hashed_nodes_at_distance_d_in_neighborhood)
        # hash the sequence of hashes of the node set at increasing
        # distances into a list of features
        hash_neighborhood = fast_hash_vec(hash_list)
        graph.nodes[root]['neigh_graph_hash'] = hash_neighborhood

    def _compute_neighborhood_graph_weight_cache(self, graph):
        assert (len(graph) > 0), 'ERROR: Empty graph'
        for u in graph.nodes():
            if graph.nodes[u].get('node', False):
                self._compute_neighborhood_graph_weight(u, graph)

    def _compute_neighborhood_graph_weight(self, root, graph):
        # list all nodes at increasing distances
        # at each distance
        # compute the arithmetic mean weight on nodes
        # compute the geometric mean weight on edges
        # compute the product of the two
        # make a list of the neighborhood_graph_weight at every distance
        neighborhood_graph_weight_list = []
        w = graph.nodes[root][self.key_weight]
        node_weight_list = np.array([w], dtype=np.float64)
        node_average = node_weight_list[0]
        edge_weight_list = np.array([1], dtype=np.float64)
        edge_average = edge_weight_list[0]
        # for all distances
        root_dist_dict = graph.nodes[root]['remote_neighbours']
        for dist in root_dist_dict.keys():
            # extract array of weights at given dist
            weight_array_at_d = np.array([graph.nodes[v][self.key_weight]
                                          for v in root_dist_dict[dist]],
                                         dtype=np.float64)
            if dist % 2 == 0:  # nodes
                node_weight_list = np.concatenate(
                    (node_weight_list, weight_array_at_d))
                node_average = np.mean(node_weight_list)
            else:  # edges
                edge_weight_list = np.concatenate(
                    (edge_weight_list, weight_array_at_d))
                edge_average = stats.gmean(edge_weight_list)
            weight = node_average * edge_average
            neighborhood_graph_weight_list.append(weight)
        graph.nodes[root]['neigh_graph_weight'] = \
            neighborhood_graph_weight_list

    def _single_vertex_breadth_first_visit(self, graph, root, max_depth):
        # the map associates to each distance value ( from 1:max_depth )
        # the list of ids of the vertices at that distance from the root
        dist_list = {}
        visited = set()  # use a set as we can end up exploring few nodes
        # q is the queue containing the frontier to be expanded in the BFV
        q = deque()
        q.append(root)
        # the map associates to each vertex id the distance from the root
        dist = {}
        dist[root] = 0
        visited.add(root)
        # add vertex at distance 0
        dist_list[0] = set()
        dist_list[0].add(root)
        while len(q) > 0:
            # extract the current vertex
            u = q.popleft()
            d = dist[u] + 1
            if d <= max_depth:
                # iterate over the neighbors of the current vertex
                for v in graph.neighbors(u):
                    if v not in visited:
                        # skip nesting edge-nodes
                        if not graph.nodes[v].get(self.key_nesting, False):
                            dist[v] = d
                            visited.add(v)
                            q.append(v)
                            if d in dist_list:
                                dist_list[d].add(v)
                            else:
                                dist_list[d] = set()
                                dist_list[d].add(v)
        graph.nodes[root]['remote_neighbours'] = dist_list

    def _compute_distant_neighbours(self, graph, max_depth):
        for n in graph.nodes():
            if graph.nodes[n].get('node', False):
                self._single_vertex_breadth_first_visit(graph, n, max_depth)

    def annotate(self,
                 graphs,
                 estimator=None,
                 reweight=1.0,
                 threshold=None,
                 scale=1,
                 vertex_features=False):
        """Return graphs with extra attributes: importance and features.

        Given a list of networkx graphs, if the given estimator is not None and
        is fitted, return a list of networkx graphs where each vertex has
        additional attributes with key 'importance' and 'weight'.
        The importance value of a vertex corresponds to the
        part of the score that is imputable to the neighborhood of radius r+d
        of the vertex. The weight value is the absolute value of importance.
        If vertex_features is True then each vertex has additional attributes
        with key 'features' and 'vector'.

        Parameters
        ----------
        estimator : scikit-learn estimator
            Scikit-learn predictor trained on data sampled from the same
            distribution. If None the vertex weights are set by default 1.

        reweight : float (default 1.0)
            The  coefficient used to weight the linear combination of the
            current weight and the absolute value of the score computed by the
            estimator.
            If reweight = 0 then do not update.
            If reweight = 1 then discard the current weight information and use
            only abs( score )
            If reweight = 0.5 then update with the arithmetic mean of the
            current weight information and the abs( score )

        threshold : float (default: None)
            If not None, threshold the importance value before computing
            the weight.

        scale : float (default: 1)
            Multiplicative factor to rescale all weights.

        vertex_features : bool (default false)
            Flag to compute the sparse vector encoding of all features that
            have that vertex as root. An attribute with key 'features' is
            created for each node that contains a CRS scipy sparse vector,
            and an attribute with key 'vector' is created that contains a
            python dictionary to store the key, values pairs.
        """
        self.estimator = estimator
        self.reweight = reweight
        self.threshold = threshold
        self.scale = scale
        self.vertex_features = vertex_features

        for graph in graphs:
            annotated_graph = self._annotate(graph)
            yield annotated_graph

    def _annotate(self, original_graph):
        # pre-processing phase: compute caches
        graph_dict = original_graph.graph
        graph = self._graph_preprocessing(original_graph)
        # extract per vertex feature representation
        data_matrix = self._compute_vertex_based_features(graph)
        # add or update weight and importance information
        graph = self._annotate_importance(graph, data_matrix)
        # add or update label information
        if self.vertex_features:
            graph = self._annotate_vector(graph, data_matrix)
        annotated_graph = _revert_edge_to_vertex_transform(graph)
        annotated_graph.graph = graph_dict
        return annotated_graph

    def _annotate_vector(self, graph, data_matrix):
        # annotate graph structure with vertex importance
        vertex_id = 0
        for v in graph.nodes():
            if graph.nodes[v].get('node', False):
                # annotate 'vector' information
                row = data_matrix.getrow(vertex_id)
                graph.nodes[v]['features'] = row
                vec_dict = {str(index): value
                            for index, value in zip(row.indices, row.data)}
                graph.nodes[v]['vector'] = vec_dict
                vertex_id += 1
        return graph

    def _compute_predictions_and_margins(self, graph, data_matrix):
        def _null_estimator_case(graph, data_matrix):
            # if we do not provide an estimator then consider default margin of
            # 1/float(len(graph)) for all vertices
            n = len(graph)
            m = data_matrix.shape[0]
            importances = np.array([1 / n] * m, dtype=np.float64)
            predictions = np.array([1] * m)
            return predictions, importances

        def _binary_classifier_case(graph, data_matrix):
            predictions = self.estimator.predict(data_matrix)
            importances = self.estimator.decision_function(data_matrix)
            return predictions, importances

        def _regression_case(graph, data_matrix):
            predicted_score = self.estimator.predict(data_matrix)
            p = [1 if v >= 0 else -1 for v in predicted_score]
            predictions = np.array(p)
            return predictions, predicted_score

        def _multiclass_case(graph, data_matrix):
            # when prediction is multiclass, use as importance the max
            # prediction
            n = len(graph)
            predictions = self.estimator.predict(data_matrix)
            predicted_score = self.estimator.decision_function(data_matrix)
            ids = np.argmax(predicted_score, axis=1)
            s = [row[col] for row, col in zip(predicted_score, ids)]
            scores = np.array(s)
            b = [self.estimator.intercept_[id] for id in ids]
            intercepts = np.array(b)
            importances = scores - intercepts + intercepts / n
            return predictions, importances

        if self.estimator is None:
            return _null_estimator_case(graph, data_matrix)
        if self.estimator.__class__.__name__ in ['SGDRegressor']:
            return _regression_case(graph, data_matrix)
        else:
            data_dim = self.estimator.intercept_.shape[0]
            if data_dim > 1:
                return _multiclass_case(graph, data_matrix)
            else:
                return _binary_classifier_case(graph, data_matrix)

    def _annotate_importance(self, graph, data_matrix):
        predictions, margins = self._compute_predictions_and_margins(
            graph, data_matrix)
        if self.threshold is not None:
            margins[margins < self.threshold] = self.threshold
        # annotate graph structure with vertex importance
        vertex_id = 0
        for v in graph.nodes():
            if graph.nodes[v].get('node', False):
                graph.nodes[v][self.key_class] = predictions[vertex_id]
                # annotate the 'importance' attribute with the margin
                graph.nodes[v][self.key_importance] = margins[vertex_id]
                # update the self.key_weight information as a linear
                # combination of the previous weight and the absolute margin
                if self.key_weight in graph.nodes[v] and self.reweight != 0:
                    graph.nodes[v][self.key_weight] = self.reweight * \
                        abs(margins[vertex_id]) +\
                        (1 - self.reweight) * \
                        graph.nodes[v][self.key_weight]
                # in case the original graph was not weighted then instantiate
                # the self.key_weight with the absolute margin
                else:
                    graph.nodes[v][self.key_weight] = abs(margins[vertex_id])
                # in all cases, rescale the weight by the scale factor
                graph.nodes[v][self.key_weight] *= self.scale
                vertex_id += 1
            if graph.nodes[v].get('edge', False):  # keep the weight of edges
                # ..unless they were unweighted, in this case add unit weight
                if self.key_weight not in graph.nodes[v]:
                    graph.nodes[v][self.key_weight] = 1
        return graph

    def _compute_vertex_based_features(self, graph):
        feature_rows = []
        for v in graph.nodes():
            # only for vertices of type 'node', i.e. not for the 'edge' type
            if graph.nodes[v].get('node', False):
                feature_list = defaultdict(lambda: defaultdict(float))
                self._transform_vertex(graph, v, feature_list)
                feature_rows.append(self._normalization(feature_list))
        data_matrix = self._convert_dict_to_sparse_matrix(feature_rows)
        return data_matrix


# -------------------------------------------------------------------

def _label_preprocessing(graph,
                         key_label='label',
                         bitmask=2 ** 20 - 1):
    for n in graph.nodes():
        label = graph.nodes[n][key_label]
        graph.nodes[n]['hlabel'] = int(hash(label) & bitmask) + 1


def _edge_to_vertex_transform(original_graph):
    """Convert edges to nodes."""
    # if operating on graphs that have already been subject to the
    # edge_to_vertex transformation, then do not repeat the transformation
    # but simply return the graph
    if 'expanded' in original_graph.graph:
        return original_graph
    else:
        # build a graph that has as vertices the original vertex set
        graph = nx.create_empty_copy(original_graph)
        graph.graph['expanded'] = True
        # add the node attribute 'node' and set it to True
        for n in graph:
            graph.nodes[n]['node'] = True
        # add a new vertex for each edge
        w = graph.number_of_nodes()
        for u, v in original_graph.edges():
            if u != v:
                graph.add_node(w)
                graph.nodes[w]['edge'] = True
                # copy the attributes of the edge
                graph.nodes[w].update(original_graph.edges[u, v])
                # and connect it to u and v
                graph.add_edge(w, u, label=None)
                graph.add_edge(w, v, label=None)
                w += 1
        return graph


def _revert_edge_to_vertex_transform(original_graph):
    """Convert nodes of type 'edge' to edges."""
    if 'expanded' in original_graph.graph:
        # start from a copy of the original graph
        graph = nx.Graph(original_graph)
        _clean_graph(graph)
        # re-wire the endpoints of edge-vertices
        for n in original_graph.nodes():
            if original_graph.nodes[n].get('edge', False):
                # extract the endpoints
                endpoints = [u for u in original_graph.neighbors(n)]
                if len(endpoints) != 2:
                    raise Exception('ERROR: more than 2 endpoints in \
                        a single edge: %s' % endpoints)
                u = endpoints[0]
                v = endpoints[1]
                # add the corresponding edge
                graph.add_edge(u, v)
                graph.edges[u, v].update(original_graph.nodes[n])
                # remove the edge-vertex
                graph.remove_node(n)
        return graph
    else:
        return original_graph


def _clean_graph(graph):
    graph.graph.pop('expanded', None)
    for n in graph.nodes():
        if graph.nodes[n].get('node', False):
            # remove stale information
            graph.nodes[n].pop('remote_neighbours', None)
            graph.nodes[n].pop('neigh_graph_hash', None)
            graph.nodes[n].pop('neigh_graph_weight', None)
            graph.nodes[n].pop('hlabel', None)


def timeit(method):
    """Time decorator."""
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        logger.debug('%s  %2.2f sec' % (method.__name__, te - ts))
        return result
    return timed


def pmap(func, iterable, chunk_size=1):
    """Multi-core map."""
    pool = mp.Pool()
    result = pool.map(func, iterable, chunksize=chunk_size)
    pool.close()
    pool.join()
    return list(result)


def ppipe(iterable, func, chunk_size=1):
    """Multi-core pipe."""
    out = pmap(func, iterable, chunk_size)
    return list(concat(out))


def configure_logging(logger, verbosity=0, filename=None):
    """Utility to configure the logging aspects.

    If filename is None then no info is stored in files.
    If filename is not None then everything that is logged is dumped to file
    (including program traces).
    Verbosity is an int that can take values: 0 -> warning,
    1 -> info, >=2 -> debug.
    All levels are displayed on stdout, not on stderr.
    Please use exceptions and asserts to output on stderr.
    """
    logger.propagate = False
    logger.handlers = []
    log_level = logging.WARNING
    if verbosity == 1:
        log_level = logging.INFO
    elif verbosity == 2:
        log_level = logging.DEBUG
    else:
        log_level = 4
    logger.setLevel(log_level)
    # create console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level)
    # create formatter
    cformatter = logging.Formatter('%(message)s')
    # add formatter to ch
    ch.setFormatter(cformatter)
    # add handlers to logger
    logger.addHandler(ch)

    if filename is not None:
        # create a file handler
        fh = logging.handlers.RotatingFileHandler(filename=filename,
                                                  maxBytes=10000000,
                                                  backupCount=10)
        fh.setLevel(logging.DEBUG)
        # create formatter
        fformatter = logging.Formatter('%(asctime)s | %(levelname)-6s | %(name)10s | %(filename)10s |\
   %(lineno)4s | %(message)s')
        # add formatter to fh
        fh.setFormatter(fformatter)
        # add handlers to logger
        logger.addHandler(fh)


def _serialize_list(items, separator='_'):
    if isinstance(items, str):
        return items
    if is_iterable(items):
        if isinstance(items, list):
            return str(separator.join([str(item) for item in items]))
        if isinstance(items, dict):
            return str(separator.join([str(key) + ':' + str(items[key])
                                       for key in items]))
    else:
        return str(items)


def serialize_dict(the_dict, full=True, offset='small'):
    """serialize_dict."""
    if the_dict:
        text = []
        for key in sorted(the_dict):
            if offset == 'small':
                line = '%10s: %s' % (key, the_dict[key])
            elif offset == 'large':
                line = '%25s: %s' % (key, the_dict[key])
            elif offset == 'very_large':
                line = '%50s: %s' % (key, the_dict[key])
            else:
                raise Exception('unrecognized option: %s' % offset)
            line = line.replace('\n', ' ')
            if full is False:
                if len(line) > 100:
                    line = line[:100] + '  ...  ' + line[-20:]
            text.append(line)
        return '\n'.join(text)
    else:
        return ""


def read(uri):
    """Abstract read function.

    EDeN can accept a URL, a file path and a python list.
    In all cases an iterable object should be returned.
    """
    if isinstance(uri, list):
        # test if it is iterable: works for lists and generators, but not for
        # strings
        return uri
    else:
        try:
                    # try if it is a URL and if we can open it
            f = requests.get(uri).text.split('\n')
        except ValueError:
            # assume it is a file object
            f = open(uri)
        return f


def is_iterable(test):
    """is_iterable."""
    if hasattr(test, '__iter__'):
        return True
    else:
        return False


def describe(data_matrix):
    """Get the shape of a sparse matrix and its average nnz."""
    return 'Instances: %3d ; Features: %d with an avg of %d per instance' % \
        (data_matrix.shape[0], data_matrix.shape[1],
         data_matrix.getnnz() / data_matrix.shape[0])


def iterator_size(iterable):
    """Length of an iterator.

    Note: if the iterable is a generator it consumes it.
    """
    if hasattr(iterable, '__len__'):
        return len(iterable)

    d = deque(enumerate(iterable, 1), maxlen=1)
    if d:
        return d[0][0]
    else:
        return 0


def random_bipartition(int_range, relative_size=.7, random_state=None):
    """random_bipartition."""
    if not random_state:
        random_state = random.random()
    random.seed(random_state)
    ids = list(range(int_range))
    random.shuffle(ids)
    split_point = int(int_range * relative_size)
    return ids[:split_point], ids[split_point:]


def selection_iterator(iterable, ids):
    """selection_iterator.

    Given an iterable and a list of ids (zero based) yield only the
    items whose id matches.
    """
    ids = sorted(ids)
    counter = 0
    for id, item in enumerate(iterable):
        if id == ids[counter]:
            yield item
            counter += 1
            if counter == len(ids):
                break


def random_bipartition_iter(iterable, relative_size=.5, random_state=1):
    """random_bipartition_iter."""
    size_iterable, iterable1, iterable2 = tee(iterable, 3)
    size = iterator_size(size_iterable)
    part1_ids, part2_ids = random_bipartition(
        size, relative_size=relative_size, random_state=random_state)
    part1_iterable = selection_iterator(iterable1, part1_ids)
    part2_iterable = selection_iterator(iterable2, part2_ids)
    return part1_iterable, part2_iterable


def store_matrix(matrix='',
                 output_dir_path='',
                 out_file_name='',
                 output_format=''):
    """store_matrix."""
    if not os.path.exists(output_dir_path):
        os.mkdir(output_dir_path)
    full_out_file_name = os.path.join(output_dir_path, out_file_name)
    if output_format == "MatrixMarket":
        if len(matrix.shape) == 1:
            raise Exception(
                "'MatrixMarket' format supports only 2D dimensional array\
                and not vectors")
        else:
            io.mmwrite(full_out_file_name, matrix, precision=None)
    elif output_format == "numpy":
        np.save(full_out_file_name, matrix)
    elif output_format == "joblib":
        joblib.dump(matrix, full_out_file_name)
    elif output_format == "text":
        with open(full_out_file_name, "w") as f:
            if len(matrix.shape) == 1:
                for x in matrix:
                    f.write("%s\n" % (x))
            else:
                raise Exception(
                    "'text' format supports only mono dimensional array\
                    and not matrices")
    logger.info("Written file: %s" % full_out_file_name)


def dump(obj, output_dir_path='', out_file_name=''):
    """dump."""
    if not os.path.exists(output_dir_path):
        os.mkdir(output_dir_path)
    full_out_file_name = os.path.join(output_dir_path, out_file_name) + ".pkl"
    joblib.dump(obj, full_out_file_name)


def load(output_dir_path='', out_file_name=''):
    """load."""
    full_out_file_name = os.path.join(output_dir_path, out_file_name) + ".pkl"
    obj = joblib.load(full_out_file_name)
    return obj


def report_base_statistics(vec, separator='\n'):
    """report_base_statistics."""
    from collections import Counter
    c = Counter(vec)
    msg = ''
    for k in c:
        msg += "class: %s count:%d (%0.2f)%s" % (
            k, c[k], c[k] / float(len(vec)), separator)
    return msg


def save_output(text=None, output_dir_path=None, out_file_name=None):
    """save_output."""
    if not os.path.exists(output_dir_path):
        os.mkdir(output_dir_path)
    full_out_file_name = os.path.join(output_dir_path, out_file_name)
    with open(full_out_file_name, 'w') as f:
        for line in text:
            f.write("%s\n" % str(line).strip())
    logger.info("Written file: %s (%d lines)" %
                (full_out_file_name, len(text)))


def sequence_dotbracket_to_graph(seq_info=None, seq_struct=None):
    """Given a sequence and the dotbracket sequence make a graph.
    Parameters
    ----------
    seq_info string
        node labels eg a sequence string
    seq_struct  string
        dotbracket string
    Returns
    -------
        returns a nx.Graph
        secondary struct associated with seq_struct
    """
    graph = nx.Graph()

    lifo = collections.defaultdict(list)
    open_brace_string={")":"(",
                "]":"[",
                ">":"<"}

    for i, (c, b) in enumerate(zip(seq_info, seq_struct)):
        graph.add_node(i, label=c, position=i)
        if i > 0:
            graph.add_edge(i, i - 1, label='-', type='backbone', len=1)
        if b in ['(','[','<']:
            lifo[b].append(i)
        if b in [')',']','>']:
            j = lifo[open_brace_string[b]].pop()
            graph.add_edge(i, j, label='=', type='basepair', len=1)

    return graph


def seq_to_graph(header, sequence):
    """Fold a sequence in a path graph."""
    seq_struct = '.' * len(sequence)
    graph = sequence_dotbracket_to_graph(seq_info=sequence,
                                         seq_struct=seq_struct)
    graph.graph['info'] = 'sequence'
    graph.graph['sequence'] = sequence
    graph.graph['structure'] = seq_struct
    graph.graph['id'] = header
    return graph


def fold(seqs):
    """Fold a list of sequences into path graphs."""
    for header, seq in seqs:
        yield seq_to_graph(header, seq)



def null_modifier(header=None, seq=None, **options):
    """Null modifier."""
    yield header, seq


def load(input, **options):
    """Load sequences."""
    return fasta_to_fasta(input, **options)

def fasta_to_fasta(input, modifier=null_modifier, **options):
    """Take a FASTA file and yield a normalised FASTA file.

    Parameters
    ----------
    input : string
        A pointer to the data source.

    normalize : bool
        If True all characters are uppercased and Ts are replaced by Us
    """
    normalize = options.get('normalize', True)
    iterable = _fasta_to_fasta(input)
    for line in iterable:
        header = line
        seq = iterable.__next__()
        if normalize:
            seq = seq.upper()
            seq = seq.replace('T', 'U')
        seqs = modifier(header=header, seq=seq, **options)
        for seq in seqs:
            yield seq


def _fasta_to_fasta(input):
    seq = ""
    with open(input, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line:
                if line[0] == '>':
                    line = line[1:]
                    if seq:
                        yield seq
                        seq = ""
                    line_str = str(line)
                    yield line_str.strip()
                else:
                    line_str = line.split()
                    if line_str:
                        seq += str(line_str[0]).strip()
    if seq:
        yield seq

