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
import networkx as nx


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


class FS_Network():

    def __init__(self, name, edges):
        """creates a network from a list of edges"""
        self.__network = nx.Graph()
        self.__name = name
        for edge in edges:
            self.__network.add_edge(edge[0], edge[1], {'score': edge[2]})

    def getGraph(self):
        return self.__network

    def addEdges(self, source, target, score):
        self.__network.add_edge(source, target, {'score':score})

    def name(self):
        """returns the name of the network"""
        return self.__name

    def edges(self):
        """returns the list of edges"""
        return self.__network.edges()

    def num_edges(self):
        """returns the number of edges in this graph"""
        return len(self.__network.edges())

    def total_score(self):
        """returns the sum of edge scores"""
        total = 0.0
        for edge in self.__network.edges(data = True):
            total += edge[2]['score']
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

    def __repr__(self):
        return "Network: %s\n# edges: %d\n" % (self.__name,
                                               len(self.__network.edges()))
    def TranslateNetwork(self, Thesaurus):
        tempG = nx.Graph()
        logging.debug("\x1b[31mNetwork:\t\x1b[0mI currently have a network of %s nodes and %s edges" % (len(self.__network.nodes()), len(self.__network.edges())))
        logging.debug("\x1b[31mNetwork:\t\x1b[0mtranslating network")
        for edge in self.__network.edges_iter(data=True):
            try:
                newS = Thesaurus[edge[0]]
                newT = Thesaurus[edge[1]]
                infoD = edge[2]
                tempG.add_edge(newS, newT, infoD)
            except KeyError:
                pass
        logging.debug("\x1b[31mNetwork:\t\x1b[0mnetwork translated...")
        logging.debug("\x1b[31mNetwork:\t\x1b[0m...into a network of %s nodes and %s edges" % (len(tempG.nodes()), len(tempG.edges())))
        
        self.__network = tempG


    def ShrinkNetwork(self, genes_all):
        tempG = self.__network.subgraph(genes_all)
        logging.debug("\x1b[31mNetwork:\t\x1b[0mnetwork shrunk...")
        logging.debug("\x1b[31mNetwork:\t\x1b[0m...into a network of %s nodes and %s edges" % (len(tempG.nodes()), len(tempG.edges())))
        
        self.__network = tempG
    


COMPUTE_NETWORK = None
ALL_GENES = None


def compute_network_scores(genes):
    """Generic method to compute network scores"""
    #network, genes, all_genes = args
    global COMPUTE_NETWORK, ALL_GENES
    network = COMPUTE_NETWORK
    all_genes = genes
    # get all edges within subnetwork
    # for edge in  G.edges.iter(all_genes, data = True)
    # 
    # if edge[0] in all_genes and edge[1] in all_genes:
    # source = edge[0]
    # target = edge[1]
    # edge[3]['Score']
    # if target not in gene_scores:
    #     gene_scores[target] = []
    # gene_scores[target].append(score)
    # if source not in gene_scores:
    #     gene_scores[source] = []
    # gene_scores[source].append(score)
    
    # also consider using H = G.subgraph(all_genes)
    # gene_scores[gene] = sum(x for x in H.edges(gene, data=True)[3]['Score'])
    
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
        final_gene_scores[gene] = -np.log(final_gene_scores[gene] + 1)
    return final_gene_scores




def compute_network_scores_FS(genes):
    """Generic method to compute network scores"""
    #network, genes, all_genes = args
    global COMPUTE_NETWORK, ALL_GENES
    network = COMPUTE_NETWORK
    all_genes = genes
    
    gene_scores = {}
    H = network.getGraph().subgraph(all_genes)
    #logging.debug("\x1b[31mNetwork:\t\x1b[0mtrimmed to a NW of %s nodes and %s edges with %s genes" %(len(H.nodes()), len(H.edges()), len(all_genes)))
    for gene in all_genes:
        gene_scores[gene] = sum(x[2]['score'] for x in H.edges(gene, data=True))
    

    final_gene_scores = {}
    for gene, scores in gene_scores.items():
        final_gene_scores[gene] = scores / len(genes)
        final_gene_scores[gene] = -np.log(final_gene_scores[gene] + 1)
    return final_gene_scores


class ScoringFunction(scoring.ScoringFunctionBase):
    """Network scoring function"""

    def __init__(self, organism, membership, matrix, weight_func=None,
                 run_in_iteration=scoring.default_network_iterations,
                 config_params=None):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, weight_func,
                                             config_params)
        self.__organism = organism
        self.__run_in_iteration = run_in_iteration
        self.__networks = None

    def name(self):
        """returns the name of this function"""
        return "Network"

    def compute(self, iteration_result, ref_matrix=None):
        """compute method"""
        iteration = iteration_result['iteration']
        if self.__run_in_iteration(iteration):
            return self.__compute()
        else:
            return None

    def __compute(self):
        """compute method, iteration is the 0-based iteration number"""
        # networks are cached
        if self.__networks == None:
            self.__networks = retrieve_networks(self.__organism)

        weight = 0.5  # TODO: for now it's fixed, we need to make them flexible
        matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                               self.gene_names())
        #network_iteration_scores = self.__create_network_iteration_scores()
        #score_means = {}  # a dictionary indexed with network names

        for network in self.__networks:
            logging.info("\x1b[31mNetwork:\t\x1b[0mCompute scores for network '%s'", network.name())
            start_time = util.current_millis()
            network_score = self.__compute_network_cluster_scores(network)
            self.__update_score_matrix(matrix, network_score, weight)
            elapsed = util.current_millis() - start_time
            logging.info("\x1b[31mNetwork:\t\x1b[0mNETWORK '%s' SCORING TIME: %f s.",
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
            map_results = pool.map(compute_network_scores_FS, args)        # Frank Schmitz - change to the FS routine
            #map_results = pool.map(compute_network_scores, args)
            pool.close()
            for cluster in xrange(1, self.num_clusters() + 1):
                result[cluster] = map_results[cluster - 1]
        else:
            for cluster in xrange(1, self.num_clusters() + 1):
                result[cluster] = compute_network_scores_FS(
                    sorted(self.rows_for_cluster(cluster)))
                
                        # Frank Schmitz - change to the FS routine
                #result[cluster] = compute_network_scores(
                #    sorted(self.rows_for_cluster(cluster)))

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
        #logging.info("Network '%s' with %d edges", network.name(),
        #             network.num_edges())
        nw_total = network.total_score()
        if nw_total > max_score:
            max_score = nw_total

    for network in networks:
        network.normalize_scores_to(max_score)
    return networks
