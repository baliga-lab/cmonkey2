# vi: sw=4 ts=4 et:
"""network.py - cMonkey network module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import numpy as np
import logging
import util
import datamatrix as dm
import scoring
import multiprocessing as mp


class Network:
    """class to represent a network graph.
    For efficiency reasons, edges is a list of [source, target, weight]
    """

    def __init__(self, name, edges, weight):
        """creates a network from a list of edges"""
        self.__name = name
        self.__edges = edges
        self.__weight = weight

    def name(self):
        """returns the name of the network"""
        return self.__name

    def edges(self):
        """returns the list of edges"""
        return self.__edges

    def weight(self):
        """returns the scoring weight of this network"""
        return self.__weight

    def num_edges(self):
        """returns the number of edges in this graph"""
        return len(self.__edges)

    def total_score(self):
        """returns the sum of edge scores"""
        total = 0.0
        for edge in self.__edges:
            total += edge[2]
        return total

    def normalize_scores_to(self, score):
        """normalizes all edge scores so that they sum up to
        the specified score"""
        total = self.total_score()
        if score != total:
            # score_e / score_total * score == score_e * (score_total / score)
            # we use this to save a division per loop iteration
            scale = float(score) / float(total)
            for edge in self.__edges:
                edge[2] = edge[2] * scale

    def edges_with_source_in(self, nodes):
        """Returns all edges containing any of the specified nodes"""
        return [edge for edge in self.__edges if edge[0] in nodes]

    def __repr__(self):
        return "Network: %s\n# edges: %d\n" % (self.__name,
                                               len(self.__edges))

    @classmethod
    def create(cls, name, edges, weight):
        """standard Factory method"""
        added = {}
        network_edges = []

        def add_edge(edge):
            """adds an edge to network_edges"""
            key = "%s:%s" % (edge[0], edge[1])
            added[key] = True
            network_edges.append(edge)

        for edge in edges:
            key = "%s:%s" % (edge[0], edge[1])
            key_rev = "%s:%s" % (edge[1], edge[0])
            if key not in added:
                add_edge(edge)
            if key_rev not in added:
                add_edge([edge[1], edge[0], edge[2]])
        return Network(name, network_edges, weight)


COMPUTE_NETWORK = None
ALL_GENES = None


def compute_network_scores(genes):
    """Generic method to compute network scores"""
    global COMPUTE_NETWORK, ALL_GENES
    network = COMPUTE_NETWORK
    all_genes = genes

    edges = network.edges_with_source_in(genes)
    fedges = [edge for edge in edges if edge[1] in all_genes]

    gene_scores = {}
    for edge in fedges:
        if edge[1] not in gene_scores:
            gene_scores[edge[1]] = []
        gene_scores[edge[1]].append(edge[2])

    final_gene_scores = {}
    for gene, scores in gene_scores.items():
        final_gene_scores[gene] = sum(scores) / len(genes)
        final_gene_scores[gene] = -np.log(final_gene_scores[gene] + 1)
    return final_gene_scores


class ScoringFunction(scoring.ScoringFunctionBase):
    """Network scoring function. Note that even though there are several
    networks, scoring can't be generalized with the default ScoringCombiner,
    since the scores are computed through weighted addition rather than
    quantile normalization"""

    def __init__(self, organism, membership, matrix, scaling_func=None,
                 run_in_iteration=scoring.default_network_iterations,
                 config_params=None):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, scaling_func,
                                             config_params)
        self.__organism = organism
        self.__run_in_iteration = run_in_iteration
        self.__networks = None
        self.__last_computed_result = None
        self.run_log = scoring.RunLog("network")
        self.__last_score_means = {}

    def name(self):
        """returns the name of this function"""
        return "Network"

    def run_logs(self):
        return [self.run_log]

    def compute(self, iteration_result, ref_matrix=None):
        """compute method
        Note: will return None if not computed yet and the result of a previous
        scoring if the function is not supposed to actually run in this iteration
        """
        iteration = iteration_result['iteration']
        if self.__run_in_iteration(iteration):
            logging.info("RUNNING A NEW NETWORK SCORING")
            self.__last_computed_result = self.__compute(iteration_result)
        self.run_log.log(self.__run_in_iteration(iteration),
                         self.scaling(iteration))
        iteration_result['networks'] = self.__last_score_means
        return self.__last_computed_result

    def __compute(self, iteration_result):
        """compute method, iteration is the 0-based iteration number"""
        # networks are cached
        if self.__networks == None:
            self.__networks = retrieve_networks(self.__organism)

        matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                               self.gene_names())

        # a dictionary that holds the scores of each gene in a given cluster
        network_iteration_scores = self.__create_network_iteration_scores()

        # a dictionary that holds the network score means for
        # each cluster, separated for each network
        score_means = {}

        for network in self.__networks:
            logging.info("Compute scores for network '%s', WEIGHT: %f",
                         network.name(), network.weight())
            start_time = util.current_millis()
            network_score = self.__compute_network_cluster_scores(network)
            self.__update_score_matrix(matrix, network_score, network.weight())
            elapsed = util.current_millis() - start_time
            logging.info("NETWORK '%s' SCORING TIME: %f s.",
                         network.name(), (elapsed / 1000.0))
            # additional scoring information, not used for the actual clustering
            score_means[network.name()] = self.__compute_cluster_score_means(
                network_score)
            self.__update_network_iteration_scores(network_iteration_scores,
                                                   network_score, network.weight())
            iteration_scores = compute_iteration_scores(
                network_iteration_scores)

        self.__last_score_means = score_means
        return matrix - matrix.quantile(0.99)

    def __compute_network_cluster_scores(self, network):
        """computes the cluster scores for the given network"""
        global COMPUTE_NETWORK, ALL_GENES
        result = {}
        use_multiprocessing = self.config_params[
            scoring.KEY_MULTIPROCESSING]
        # Set the huge memory objects into globals
        # These are readonly anyways, but using Manager.list() or something
        # similar brings this down to a crawl
        COMPUTE_NETWORK = network
        ALL_GENES = self.gene_names()

        if use_multiprocessing:
            pool = mp.Pool()
            args = []
            for cluster in xrange(1, self.num_clusters() + 1):
                args.append(sorted(self.rows_for_cluster(cluster)))
            map_results = pool.map(compute_network_scores, args)
            pool.close()
            for cluster in xrange(1, self.num_clusters() + 1):
                result[cluster] = map_results[cluster - 1]
        else:
            for cluster in xrange(1, self.num_clusters() + 1):
                result[cluster] = compute_network_scores(
                    sorted(self.rows_for_cluster(cluster)))

        return result

    def __update_score_matrix(self, matrix, network_score, weight):
        """add values into the result score matrix"""
        for cluster in xrange(1, self.num_clusters() + 1):
            for row_index in xrange(self.matrix().num_rows()):
                gene = self.gene_at(row_index)
                if gene in network_score[cluster].keys():
                    weighted_score = network_score[cluster][gene] * weight
                    matrix[row_index][cluster - 1] += weighted_score

    # The functions below are computed by cMonkey for stats, we don't
    # use them right now, but keep them around for debugging and
    # integration of stats functionality
    def __create_network_iteration_scores(self):
        """creates initialized network iteration scores"""
        result = {}
        for cluster in xrange(1, self.num_clusters() + 1):
            result[cluster] = {}
        return result

    def __compute_cluster_score_means(self, network_score):
        """compute the score means on the given network score"""
        result = {}
        for cluster in xrange(1, self.num_clusters() + 1):
            cluster_scores = []
            for gene in sorted(self.rows_for_cluster(cluster)):
                if gene in network_score[cluster].keys():
                    cluster_scores.append(network_score[cluster][gene])
                else:
                    cluster_scores.append(0.0)
            result[cluster] = util.trim_mean(cluster_scores, 0.05)
        return result

    def __update_network_iteration_scores(self, result, network_score, weight):
        """compute network iteration scores"""
        for cluster in xrange(1, self.num_clusters() + 1):
            for gene in sorted(self.rows_for_cluster(cluster)):
                if gene not in result[cluster].keys():
                    result[cluster][gene] = 0.0
                if gene in network_score[cluster].keys():
                    weighted_score = network_score[cluster][gene] * weight
                    result[cluster][gene] += weighted_score
        return result


def compute_iteration_scores(network_iteration_scores):
    """called 'cluster.ns' in the original cMonkey"""
    result = {}
    for cluster in network_iteration_scores:
        cluster_scores = []
        for _, score in network_iteration_scores[cluster].items():
            cluster_scores.append(score)
        result[cluster] = util.trim_mean(cluster_scores, 0.05)
    return result


def retrieve_networks(organism):
    """retrieves the networks provided by the organism object and
    possibly other sources, doing some normalization if necessary
    Note: wanted to make it private, but the scoring function
    can not see it after doing so"""
    networks = organism.networks()
    max_score = 0
    for network in networks:
        #logging.info("Network '%s' with %d edges", network.name(),
        #             network.num_edges())
        nw_total = network.total_score()
        if nw_total > max_score:
            max_score = nw_total

    for network in networks:
        network.normalize_scores_to(max_score)
    return networks
