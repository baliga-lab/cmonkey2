"""network.py - cMonkey network module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import numpy
import logging
import util
import membership as memb


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


def compute_network_scores(network, genes, all_genes):
    """Generic method to compute network scores
    TODO: maybe should be part of Network class"""
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


class ScoringFunction(memb.ScoringFunctionBase):
    """Network scoring function"""

    def __init__(self, organism, membership, matrix, weight_func=None):
        """Create scoring function instance"""
        memb.ScoringFunctionBase.__init__(self, membership,
                                          matrix, weight_func)
        self.__organism = organism

    def compute(self, iteration):
        """compute method, iteration is the 0-based iteration number"""

        result = {}  # a dictionary indexed with network names
        networks = self.retrieve_networks(self.__organism)

        # TODO: here add the network scores weighted (0.5 for the two networks
        # for now)
        weight = 0.5  # for now it's fixed, we need to make them flexible
        network_iteration_scores = {}
        for cluster in range(1, self.num_clusters() + 1):
            network_iteration_scores[cluster] = {}

        for network in networks:
            logging.info("Compute scores for network '%s'", network.name())
            network_score = {}
            cluster_score_means = {}

            for cluster in range(1, self.num_clusters() + 1):
                cluster_genes = sorted(self.rows_for_cluster(cluster))
                network_score[cluster] = compute_network_scores(
                    network, cluster_genes, self.gene_names())
                # build network scoring based on cluster membership
                cluster_scores = []
                for gene in sorted(self.rows_for_cluster(cluster)):
                    # init iteration score if non-existent
                    if gene not in network_iteration_scores[cluster].keys():
                        network_iteration_scores[cluster][gene] = 0.0

                    if gene in network_score[cluster].keys():
                        cluster_scores.append(network_score[cluster][gene])
                        network_iteration_scores[cluster][gene] += network_score[cluster][gene] * weight
                    else:
                        cluster_scores.append(0.0)

                cluster_score_means[cluster] = util.trim_mean(cluster_scores,
                                                              0.05)
            result[network.name()] = cluster_score_means

        iteration_scores = {}
        for cluster in network_iteration_scores:
            cluster_scores = []
            for gene, score in network_iteration_scores[cluster].items():
                cluster_scores.append(score)
            iteration_scores[cluster] = util.trim_mean(cluster_scores, 0.05)
        result['all'] = iteration_scores
        return result

    def retrieve_networks(self, organism):
        """retrieves the networks provided by the organism object and
        possibly other sources, doing some normalization if necessary"""
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
