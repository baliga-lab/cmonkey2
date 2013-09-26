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
import cPickle
import os.path


class Network:
    """class to represent a network graph.
    The graph is considered undirected
    For efficiency reasons, edges is a list of [source, target, weight]
    """

    def __init__(self, name, edges, weight, dummy):
        """creates a network from a list of edges"""
        self.name = name
        self.edges = edges
        self.weight = weight
        self.__compute_edges_with_source()

    def __compute_edges_with_source(self):
        self.edges_with_source = {}
        for edge in self.edges:
            if edge[0] not in self.edges_with_source:
                self.edges_with_source[edge[0]] = []
            if edge[1] not in self.edges_with_source:
                self.edges_with_source[edge[1]] = []
            self.edges_with_source[edge[0]].append(edge)
            self.edges_with_source[edge[1]].append(edge)

    def validate(self, synonyms, genes):
        # remap first
        new_edges = []
        for n0, n1, score in self.edges:
            n0 = synonyms[n0] if n0 in synonyms else n0
            n1 = synonyms[n1] if n1 in synonyms else n1
            new_edges.append((n0, n1, score))
        self.edges = new_edges
        self.__compute_edges_with_source()

        # then validate
        found = []
        for g in genes:
            primary = synonyms.get(g, g)
            for n0, n1, score in self.edges:
                if primary == n0 or primary == n1:
                    found.append(primary)
        if len(found) < len(genes) / 2:
            print edges
            raise(Exception("only %d genes found in edges" % len(found)))

    def num_edges(self):
        """returns the number of edges in this graph"""
        return len(self.edges)

    def total_score(self):
        """returns the sum of edge scores"""
        return sum(edge[2] for edge in self.edges) * 2

    def normalize_scores_to(self, score):
        """normalizes all edge scores so that they sum up to
        the specified score"""
        total = self.total_score()
        if score != total:
            # score_e / score_total * score == score_e * (score_total / score)
            # we use this to save a division per loop iteration
            scale = float(score) / float(total)
            for edge in self.edges:
                edge[2] = edge[2] * scale

    def edges_with_node(self, node):
        """returns the edges where node is a node of"""
        if node in self.edges_with_source:
            return self.edges_with_source[node]
        else:
            return []

    def __repr__(self):
        return "Network: %s\n# edges: %d\n" % (self.name,
                                               len(self.edges))

    @classmethod
    def create(cls, name, edges, weight, check_size=True):
        """standard Factory method"""
        if edges == None:
            raise Exception("no edges specified in network '%s'" % name)
        added = {}
        network_edges = []

        for edge in edges:
            key = "%s:%s" % (edge[0], edge[1])
            key_rev = "%s:%s" % (edge[1], edge[0])
            if key not in added and key_rev not in added:
                network_edges.append(edge)
            added[key] = True
            added[key_rev] = True

        if check_size and len(network_edges) < 10:
            raise Exception("Error: only %d edges in network '%s'" % (len(network_edges), name))
        logging.info("Created network '%s' with %d edges", name, len(network_edges))
        return Network(name, network_edges, weight, 0)


COMPUTE_NETWORK = None
ALL_GENES = None
NETWORK_SCORE_MEMBERSHIP = None

def compute_network_scores(cluster):
    """Generic method to compute network scores"""
    global COMPUTE_NETWORK, ALL_GENES, NETWORK_SCORE_MEMBERSHIP
    network = COMPUTE_NETWORK

    genes = NETWORK_SCORE_MEMBERSHIP.rows_for_cluster(cluster)
    gene_scores = {}
    for gene in genes:
        edges = network.edges_with_node(gene)
        for edge in edges:
            other_gene = edge[0]
            if other_gene == gene:
                other_gene = edge[1]
            if other_gene in ALL_GENES:
                if other_gene not in gene_scores:
                    gene_scores[other_gene] = []
                gene_scores[other_gene].append(edge[2])


    final_gene_scores = {}
    for gene, scores in gene_scores.items():
        final_gene_scores[gene] = sum(scores) / len(genes)
        final_gene_scores[gene] = -np.log(final_gene_scores[gene] + 1)
    return final_gene_scores


def compute_mean(score_means):
    means = {}
    for network, cluster_score_means in score_means.items():
        total = 0.0
        for score in cluster_score_means.values():
            total = total + score
            means[network] = total / len(cluster_score_means)
    return means


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
                                             run_in_iteration,
                                             config_params)
        self.__organism = organism
        self.__networks = None
        self.__last_computed_result = None
        self.run_log = scoring.RunLog("network", config_params)

    def name(self):
        """returns the name of this function"""
        return "Network"

    def run_logs(self):
        return [self.run_log]

    def network_scores_pickle_path(self):
        return "%s/%s_scores_last.pkl" % (self.config_params['output_dir'],
                                          self.name())

    def current_network_scores(self):
        if self.network_scores != None:
            return self.network_scores
        else:
            with open(self.network_scores_pickle_path()) as infile:
                result = cPickle.load(infile)
            return result

    def store_checkpoint_data(self, shelf):
        """Default implementation does not store checkpoint data"""
        shelf['network.scores'] = self.network_scores

    def restore_checkpoint_data(self, shelf):
        """Default implementation does not store checkpoint data"""
        self.network_scores = shelf['network.scores']

    def compute(self, iteration_result, ref_matrix=None):
        """overridden compute for storing additional information"""
        result = scoring.ScoringFunctionBase.compute(self, iteration_result, ref_matrix)
        iteration_result['networks'] = self.__update_score_means(
            self.current_network_scores())
        self.network_scores = None
        return result

    def compute_force(self, iteration_result, ref_matrix=None):
        """overridden compute for storing additional information"""
        result = scoring.ScoringFunctionBase.compute_force(self, iteration_result, ref_matrix)
        iteration_result['networks'] = self.__update_score_means(
            self.current_network_scores())
        self.network_scores = None
        return result

    def networks(self):
        """networks are cached"""
        if self.__networks == None:
            self.__networks = retrieve_networks(self.__organism)
            if self.config_params['remap_network_nodes']:
                # network names are non-primary, this can happen
                # when the user makes up their own data
                for network in self.__networks:
                    network.validate(self.__organism.thesaurus(),
                                     self.gene_names())
        return self.__networks

    def __update_score_means(self, network_scores):
        """returns the score means, adjusted to the current cluster setup"""
        # a dictionary that holds the network score means for
        # each cluster, separated for each network
        score_means = {}
        for network in self.networks():
            score_means[network.name] = self.__compute_cluster_score_means(
                network_scores[network.name])
        return compute_mean(score_means)

    def do_compute(self, iteration_result, ref_matrix=None):
        """compute method, iteration is the 0-based iteration number"""

        matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                               self.gene_names())

        # a dictionary that holds the scores of each gene in a given cluster
        network_iteration_scores = {cluster: {}
                                    for cluster in xrange(1, self.num_clusters() + 1)}
        network_scores = {}
        for network in self.networks():
            logging.info("Compute scores for network '%s', WEIGHT: %f",
                         network.name, network.weight)
            start_time = util.current_millis()
            network_score = self.__compute_network_cluster_scores(network)
            network_scores[network.name] = network_score
            self.__update_score_matrix(matrix, network_score, network.weight)
            elapsed = util.current_millis() - start_time
            logging.info("NETWORK '%s' SCORING TIME: %f s.",
                         network.name, (elapsed / 1000.0))
            # additional scoring information, not used for the actual clustering
            self.__update_network_iteration_scores(network_iteration_scores,
                                                   network_score, network.weight)
            iteration_scores = compute_iteration_scores(network_iteration_scores)

        outfile_path = self.network_scores_pickle_path()
        outfile_dir = os.path.dirname(outfile_path)
        if not os.path.isdir(outfile_dir):
            from os import makedirs
            makedirs(outfile_dir)
        with open(self.network_scores_pickle_path(), 'w') as outfile:
            cPickle.dump(network_scores, outfile)
        # immediately use in means computation
        self.network_scores = network_scores
        matrix.subtract_with_quantile(0.99)
        return matrix

    def __compute_network_cluster_scores(self, network):
        """computes the cluster scores for the given network"""
        global COMPUTE_NETWORK, ALL_GENES, NETWORK_SCORE_MEMBERSHIP
        result = {}
        use_multiprocessing = self.config_params[
            scoring.KEY_MULTIPROCESSING]
        # Set the huge memory objects into globals
        # These are readonly anyways, but using Manager.list() or something
        # similar brings this down to a crawl
        COMPUTE_NETWORK = network
        ALL_GENES = set(self.gene_names())  # optimization: O(1) lookup
        NETWORK_SCORE_MEMBERSHIP = self.membership()

        if use_multiprocessing:
            pool = mp.Pool()
            map_results = pool.map(compute_network_scores, xrange(1, self.num_clusters() + 1))
            pool.close()
            pool.join()
            for cluster in xrange(1, self.num_clusters() + 1):
                result[cluster] = map_results[cluster - 1]
        else:
            for cluster in xrange(1, self.num_clusters() + 1):
                result[cluster] = compute_network_scores(cluster)
        # cleanup
        COMPUTE_NETWORK = None
        ALL_GENES = None
        NETWORK_SCORE_MEMBERSHIP = None
        return result

    def __update_score_matrix(self, matrix, network_score, weight):
        """add values into the result score matrix"""
        mvalues = matrix.values
        gene_names = self.gene_names()
        for cluster in xrange(1, self.num_clusters() + 1):
            cluster_genes = set(network_score[cluster].keys())
            for row_index in xrange(self.matrix().num_rows):
                gene = gene_names[row_index]
                if gene in cluster_genes:
                    weighted_score = network_score[cluster][gene] * weight
                    mvalues[row_index][cluster - 1] += weighted_score

    def __compute_cluster_score_means(self, network_score):
        """compute the score means on the given network score"""
        result = {}
        for cluster in xrange(1, self.num_clusters() + 1):
            ## TODO: we can rewrite this as a list comprehension. Check whether that
            ## will be faster
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
