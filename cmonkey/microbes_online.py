"""microbes_online.py - operon prediction module cMonkey
For organisms that have operons, cMonkey can read use operon predictions
to improve its scoring. This module defines access to prediction data through
various services.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import logging
from util import read_url_cached, DelimitedFile
from network import Network, NetworkEdge


MICROBES_ONLINE_BASE_URL = 'http://www.microbesonline.org'


class MicrobesOnline:
    """Interface to Microbes Online web service"""

    def __init__(self, base_url=MICROBES_ONLINE_BASE_URL,
                 cache_dir='cache'):
        """creates a MicrobesOnline service instance"""
        self.base_url = base_url
        self.cache_dir = cache_dir

    def get_operon_predictions_for(self, organism_id):
        """Retrieve operon predictions for the specified organism"""
        logging.info("MicrobesOnline.get_operon_predictions_for(%s)",
                     organism_id)
        url = '/'.join([self.base_url, 'operons',
                       'gnc%s.named' % str(organism_id)])
        cache_file = '/'.join([self.cache_dir,
                              'gnc%s.named' % str(organism_id)])
        return read_url_cached(url, cache_file)


def make_operon_edges(operon, organism, features):
    """take an operon as a list of gene names, determines the head out of
    these gene names and generates edges from the head to each gene in the
    operon.
    The head is is determined as follows:
    1. retrieve the gene coordinates for each gene in the operon
    2. if most genes are on the forward strand, the head is the one
       with the lowest start position
    3. if most genes are on the reverse strand, the head is the one
       with the highest end position
    This function returns an empty result if
    1. the same amount of genes are on the forward and reverse strand
    2. the gene coordinates can't be retrieved
    """
    def get_reverse_head(feature_map):
        """determine reverse head of the operon"""
        max_gene = None
        max_end = 0
        for (gene, feature) in feature_map.items():
            if feature.end() > max_end:
                max_end = feature.end()
                max_gene = gene
        return max_gene

    def get_forward_head(feature_map):
        """determine forward head of the operon"""
        min_gene = None
        min_start = sys.maxint
        for (gene, feature) in feature_map.items():
            if feature.start() < min_start:
                min_start = feature.start()
                min_gene = gene
        return min_gene

    feature_map = {}
    num_reverse = 0
    for gene in operon:
        feature_map[gene] = features[organism.feature_id_for(gene)]
        if feature_map[gene].is_reverse():
            num_reverse += 1
    num_total = len(operon)
    percent_reverse = float(num_reverse) / float(num_total)
    if percent_reverse > 0.6:
        head = get_reverse_head(feature_map)
    elif percent_reverse < 0.4:
        head = get_forward_head(feature_map)
    else:
        logging.warning("can't determine head of operon - amounts " +
                        "of reverse and forward genes are too similar (%f-%f)",
                        percent_reverse, 1.0 - percent_reverse)
        return []
    return [(head, gene) for gene in operon]


def build_operons(names1, names2):
    """build the list of operons given two name lists"""
    def first_row_containing(alist, aname):
        """returns in a list containing the specified name"""
        for row in alist:
            if aname in row:
                return row
        return None

    operons = []
    for i in range(len(names1)):
        found = first_row_containing(operons, names1[i])
        if found:
            found.append(names2[i])
        else:
            operons.append([names1[i], names2[i]])
    return operons


def make_edges_from_predictions(predictions, organism):
    """turn a list of predictions into a list of network edges"""
    def build_names():
        """builds the gene name lists from the predictions"""
        names1 = []
        names2 = []
        for prediction in predictions:
            if prediction[0] not in names1:
                names1.append(prediction[0])
            if prediction[1] not in names2:
                names2.append(prediction[1])
        return names1, names2

    names1, names2 = build_names()
    features = organism.features_for_genes(names1 + names2)
    operons = build_operons(names1, names2)
    edges = []
    for operon in operons:
        edges.extend(make_operon_edges(operon, organism, features))
    return edges


def get_network_factory(microbes_online):
    """function to create a network factory method"""
    def make_network(organism):
        """factory method to create a network from operon predictions"""
        logging.info("MicrobesOnline - make_network()")
        preds_text = microbes_online.get_operon_predictions_for(
            organism.taxonomy_id())
        dfile = DelimitedFile.create_from_text(preds_text, has_header=True)
        preds = [(line[2], line[3]) for line in dfile.lines()
                 if line[6] == 'TRUE']
        pred_edges = make_edges_from_predictions(preds, organism)
        edges = [NetworkEdge(edge[0], edge[1], 1000) for edge in pred_edges]
        return Network.create('operons', edges)

    return make_network


__all__ = ['MicrobesOnline', 'get_network_factory']
