'''
Created on Apr 2, 2012

@author: frank
'''

import numpy as np
import logging
import util
import datamatrix as dm
import scoring
import multiprocessing as mp
import networkx as nx
import network as nwGeneric

# TODO:
# add networks for
# GO
# pathways

COMPUTE_NETWORK = None
ALL_GENES = None



class Network(nwGeneric.Network):

    def __init__(self, name, weight):
        """creates a network from a list of edges"""
        self.__network = None
        self.__name = name
        self.__weight = weight
    def buildNWfromList(self, edges):
        if self.__network == None:
            self.__network = nx.Graph()
        for edge in edges:
            self.__network.add_edge(edge[0], edge[1], {'score': edge[2]})
        logging.debug("\x1b[31mNetwork:\t\x1b[0mbuilt network from edgelist")

    def getGraph(self):
        return self.__network
    def putGraph(self, graph):
        self.__network = graph
        logging.debug("\x1b[31mNetwork:\t\x1b[0mimported network from external graph")
    def addEdges(self, source, target, score):
        self.__network.add_edge(source, target, {'score':score})
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
    






class ScoringFunction(nwGeneric.ScoringFunction):
    """Network scoring function. Note that even though there are several
    networks, scoring can't be generalized with the default ScoringCombiner,
    since the scores are computed through weighted addition rather than
    quantile normalization"""
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
            map_results = pool.map(compute_network_scores, args)     # Frank Schmitz
            pool.close()
            for cluster in xrange(1, self.num_clusters() + 1):
                result[cluster] = map_results[cluster - 1]
        else:
            for cluster in xrange(1, self.num_clusters() + 1):
                result[cluster] = compute_network_scores(            # Frank Schmitz
                    sorted(self.rows_for_cluster(cluster)))
        return result



def compute_network_scores(genes):
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


