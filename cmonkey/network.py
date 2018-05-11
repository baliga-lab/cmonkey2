# vi: sw=4 ts=4 et:
"""network.py - cMonkey network module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import numpy as np
import logging
import os.path

import cmonkey.util as util
import cmonkey.datamatrix as dm
import cmonkey.scoring as scoring

# Python2/Python3 compatibility
try:
    xrange
except NameError:
    xrange = range


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
        """Change the names in the network to have the standard names in the
            synonyms (elswhere call the thesaurus).  Problem: it does not
            also rename the ratios matrix to the standard names

             Keyword arguments:
             synonyms  -- The thesaurus.
             genes     -- The gene names from the ratios.

             Usage:
             self.validate(synonyms, genes)
        """
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
            print(edges)
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
            self.edges = [(edge[0], edge[1], edge[2] * scale) for edge in self.edges]
        self.__compute_edges_with_source()

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
    def create(cls, name, edges, weight, organism=None, ratios=None,
               check_size=True):
        """standard Factory method"""
        logging.debug("Network.create() called with %d edges", len(edges))
        if edges is None:
            raise Exception("no edges specified in network '%s'" % name)
        added = set([])
        network_edges = []
        nodes = set()
        for edge in edges:
            nodes.add(edge[0])
            nodes.add(edge[1])
        """Shrink the number of edges to the ones that are actually usable. These
        are selected by the following considerations:
        # 1. check nodes that are in the thesaurus
        # 2. check gene names that are in the ratios matrix, but not in the network
        # 3. keep the nodes that are in the ratios and are in the thesaurus
        """
        num_nodes_orig = len(nodes)
        if organism:
            thesaurus = organism.thesaurus()
            nodes = {n for n in nodes if n in thesaurus}
            if ratios:
                cano_nodes = {thesaurus[n] for n in nodes}
                cano_genes = {thesaurus[row] for row in ratios.row_names
                              if row in thesaurus}
                probes_in = [gene for gene in cano_genes if gene in cano_nodes]
                nodes = {n for n in nodes if thesaurus[n] in probes_in}

        logging.debug("# nodes in network '%s': %d (of %d)", name, len(nodes), num_nodes_orig)

        for edge in edges:
            # we ignore self-edges, and edges with nodes not in the final nodes
            if edge[0] != edge[1] and edge[0] in nodes and edge[1] in nodes:
                key = "%s:%s" % (edge[0], edge[1])
                key_rev = "%s:%s" % (edge[1], edge[0])
                if key not in added and key_rev not in added:
                    network_edges.append((edge[0], edge[1], edge[2]))
                added.add(key)
                added.add(key_rev)

        if check_size and len(network_edges) < 10:
            raise Exception("Error: only %d edges in network '%s'" % (len(network_edges), name))
        logging.debug("Created network '%s' with %d edges", name, len(network_edges))
        return Network(name, network_edges, weight, 0)


COMPUTE_NETWORK = None
ALL_GENES = None
NETWORK_SCORE_MEMBERSHIP = None


def compute_network_scores(cluster):
    """Generic method to compute network scores"""
    global COMPUTE_NETWORK, ALL_GENES, NETWORK_SCORE_MEMBERSHIP
    network = COMPUTE_NETWORK

    genes = sorted(NETWORK_SCORE_MEMBERSHIP.rows_for_cluster(cluster))
    gene_scores = {}

    for gene in genes:
        # TODO: optimization: we can use numpy arrays for the scores array
        # and then sum
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



class ScoringFunction(scoring.ScoringFunctionBase):
    """Network scoring function. Note that even though there are several
    networks, scoring can't be generalized with the default ScoringCombiner,
    since the scores are computed through weighted addition rather than
    quantile normalization"""

    def __init__(self, function_id, cmrun):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, function_id, cmrun)
        self.__networks = None
        self.run_log = scoring.RunLog(function_id, cmrun.dbsession(),
                                      self.config_params)

    def initialize(self, args):
        """process additional parameters"""
        self.weights = {nw['type']: nw['weight'] for nw in args['networks']}

    def run_logs(self):
        return [self.run_log]

    def compute(self, iteration_result, ref_matrix=None):
        """overridden compute for storing additional information"""
        result = scoring.ScoringFunctionBase.compute(self, iteration_result, ref_matrix)
        iteration_result['networks'] = self.score_means

        return result

    def compute_force(self, iteration_result, ref_matrix=None):
        """overridden compute for storing additional information"""
        result = scoring.ScoringFunctionBase.compute_force(self, iteration_result, ref_matrix)
        iteration_result['networks'] = self.score_means

        return result

    def networks(self):
        """networks are cached"""
        if self.__networks is None:
            self.__networks = retrieve_networks(self.organism)
            if self.config_params['remap_network_nodes']:
                # network names are non-primary, this can happen
                # when the user makes up their own data
                for network in self.__networks:
                    network.validate(self.organism.thesaurus(),
                                     self.gene_names())
        return self.__networks

    def __update_score_means(self, network_scores):
        """returns the score means, adjusted to the current cluster setup"""
        # a dictionary that holds the network score means for
        # each cluster, separated for each network
        if network_scores:
            score_means = {network.name: self.__compute_cluster_score_means(network_scores[network.name])
                           for network in self.networks()}
            return {network: np.average(np.array(list(cluster_score_means.values())))
                    for network, cluster_score_means in score_means.items()}
        return {}

    def do_compute(self, iteration_result, ref_matrix=None):
        """compute method, iteration is the 0-based iteration number"""

        matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                               self.gene_names())
        network_scores = {}
        for network in self.networks():
            logging.debug("Compute scores for network '%s', WEIGHT: %f",
                          network.name, network.weight)
            start_time = util.current_millis()
            network_score = self.__compute_network_cluster_scores(network)
            network_scores[network.name] = network_score
            self.__update_score_matrix(matrix, network_score, network.weight)
            elapsed = util.current_millis() - start_time
            logging.debug("NETWORK '%s' SCORING TIME: %f s.",
                          network.name, (elapsed / 1000.0))

        # compute and store score means
        self.score_means = self.__update_score_means(network_scores)
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
        NETWORK_SCORE_MEMBERSHIP = self.membership

        if use_multiprocessing:
            with util.get_mp_pool(self.config_params) as pool:
                map_results = pool.map(compute_network_scores, xrange(1, self.num_clusters() + 1))
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
            for row_index in xrange(self.ratios.num_rows):
                gene = gene_names[row_index]
                if gene in cluster_genes:
                    weighted_score = network_score[cluster][gene] * weight
                    mvalues[row_index][cluster - 1] += weighted_score

    def __compute_cluster_score_means(self, network_score):
        """compute the score means on the given network score"""
        result = {}
        for cluster in xrange(1, self.num_clusters() + 1):
            cluster_scores = [network_score[cluster][gene]
                              if gene in network_score[cluster] else 0.0
                              for gene in self.rows_for_cluster(cluster)]
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
        #logging.debug("Network '%s' with %d edges", network.name(),
        #              network.num_edges())
        nw_total = network.total_score()
        if nw_total > max_score:
            max_score = nw_total

    for network in networks:
        network.normalize_scores_to(max_score)
    return networks
