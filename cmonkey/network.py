"""network.py - cMonkey network module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import numpy
import logging
import util
import datamatrix as dm
import scoring
import multiprocessing as mp


class NetworkEdge:
    """class to represent a network edge"""

    def __init__(self, source, target, score):
        """create an edge instance"""
        self.__source = source
        self.__target = target
        self.__score = score

    def source(self):
        """returns the source node"""
        return self.__source

    def target(self):
        """returns the target node"""
        return self.__target

    def score(self):
        """returns the edge score"""
        return self.__score

    def set_score(self, score):
        """sets a new score for this edge"""
        self.__score = score

    def source_in(self, nodes):
        """checks whether this edge's source is in the specified nodes"""
        return self.__source in nodes

    def target_in(self, nodes):
        """checks whether this edge's target is in the specified nodes"""
        return self.__target in nodes

    def __str__(self):
        """returns string representation"""
        return "%s -> %s w = %s" % (self.__source, self.__target,
                                    str(self.__score))

    def __repr__(self):
        """returns string representation"""
        return str(self)


class Network:
    """class to represent a network graph"""

    def __init__(self, name, edges):
        """creates a network from a list of edges"""
        self.__name = name
        self.__edges = edges

    def name(self):
        """returns the name of the network"""
        return self.__name

    def edges(self):
        """returns the list of edges"""
        return self.__edges

    def num_edges(self):
        """returns the number of edges in this graph"""
        return len(self.__edges)

    def total_score(self):
        """returns the sum of edge scores"""
        total = 0.0
        for edge in self.__edges:
            total += edge.score()
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
                edge.set_score(edge.score() * scale)

    def edges_with_source_in(self, nodes):
        """Returns all edges containing any of the specified nodes"""
        return [edge for edge in self.__edges if edge.source_in(nodes)]

    def __repr__(self):
        return "Network: %s\n# edges: %d\n" % (self.__name,
                                               len(self.__edges))

    @classmethod
    def create(cls, name, edges):
        """standard Factory method"""
        added = {}
        network_edges = []

        def add_edge(edge):
            """adds an edge to network_edges"""
            key = "%s:%s" % (edge.source(), edge.target())
            added[key] = True
            network_edges.append(edge)

        for edge in edges:
            key = "%s:%s" % (edge.source(), edge.target())
            key_rev = "%s:%s" % (edge.target(), edge.source())
            if key not in added:
                add_edge(edge)
            if key_rev not in added:
                add_edge(NetworkEdge(edge.target(), edge.source(),
                                     edge.score()))
        return Network(name, network_edges)


def compute_network_scores(args):
    """Generic method to compute network scores"""
    network, genes, all_genes = args
    edges = network.edges_with_source_in(genes)
    fedges = [edge for edge in edges if edge.target_in(all_genes)]

    gene_scores = {}
    for edge in fedges:
        if edge.target() not in gene_scores:
            gene_scores[edge.target()] = []
        gene_scores[edge.target()].append(edge.score())

    final_gene_scores = {}
    for gene, scores in gene_scores.items():
        final_gene_scores[gene] = sum(scores) / len(genes)
        final_gene_scores[gene] = -numpy.log(final_gene_scores[gene] + 1)
    return final_gene_scores


class ScoringFunction(scoring.ScoringFunctionBase):
    """Network scoring function"""

    def __init__(self, organism, membership, matrix, weight_func=None,
                 interval=0, config_params=None):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, weight_func,
                                             config_params)
        self.__organism = organism
        self.__interval = interval

    def name(self):
        """returns the name of this function"""
        return "Network"

    def compute(self, iteration, ref_matrix=None):
        """compute method"""

        if (self.__interval == 0 or
            (iteration > 0 and (iteration % self.__interval == 0))):
            return self.__compute()
        else:
            return None

    def __compute(self):
        """compute method, iteration is the 0-based iteration number"""
        networks = retrieve_networks(self.__organism)
        weight = 0.5  # for now it's fixed, we need to make them flexible
        matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                               self.gene_names())
        #network_iteration_scores = self.__create_network_iteration_scores()
        #score_means = {}  # a dictionary indexed with network names

        for network in networks:
            logging.info("Compute scores for network '%s'", network.name())
            start_time = util.current_millis()
            network_score = self.__compute_network_cluster_scores(network)
            self.__update_score_matrix(matrix, network_score, weight)
            elapsed = util.current_millis() - start_time
            logging.info("NETWORK '%s' SCORING TIME: %f s.",
                         network.name(), (elapsed / 1000.0))
            #score_means[network.name()] = self.__compute_cluster_score_means(
            #    network_score)
            #self.__update_network_iteration_scores(network_iteration_scores,
            #                                       network_score, weight)
            #iteration_scores = __compute_iteration_scores(
            #    network_iteration_scores)
        return matrix - matrix.quantile(0.99)

    def __compute_network_cluster_scores(self, network):
        """computes the cluster scores for the given network"""
        result = {}
        use_multiprocessing = self.config_params[
            scoring.KEY_MULTIPROCESSING]
        if use_multiprocessing:
            pool = mp.Pool()
            args = []
            for cluster in xrange(1, self.num_clusters() + 1):
                args.append((network, sorted(self.rows_for_cluster(cluster)),
                    self.gene_names()))
            map_results = pool.map(compute_network_scores, args)
            pool.close()
            for cluster in xrange(1, self.num_clusters() + 1):
                result[cluster] = map_results[cluster - 1]
        else:
            for cluster in xrange(1, self.num_clusters() + 1):
                result[cluster] = compute_network_scores(
                    (network, sorted(self.rows_for_cluster(cluster)),
                    self.gene_names()))

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


def __compute_iteration_scores(network_iteration_scores):
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
        logging.info("Network '%s' with %d edges", network.name(),
                     network.num_edges())
        nw_total = network.total_score()
        if nw_total > max_score:
            max_score = nw_total

    for network in networks:
        network.normalize_scores_to(max_score)
    return networks
