"""stringdb.py - cMonkey STRING database interface module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import logging
import re
import math
import util
import network
import patches


STRING_FILE2 = 'string_links_64091.tab'
PROTEIN_PREFIX = re.compile('^string:\d+[.]')

def normalize_edges_to_max_score(edges, max_score):
    """normalize scores to 1000, for combined scores"""
    def normalize(edge_score):
        score = edge_score / max_score * 1000.0
        return 1000 * math.exp(score / 1000.0) / math.exp(1.0)

    return [(edge[0], edge[1], normalize(edge[2])) for edge in edges]


def get_network_factory2(organism_code, filename, weight, sep='\t',
                         normalized=False):
    """STRING network factory from preprocessed edge file
    (protein1, protein2, combined_score), scores are already
    normalized to 1000.
    This is the standard factory method used for Microbes.
    """
    def read_edges2(filename):
        """just read a preprocessed file, much faster to debug"""
        logging.info("stringdb.read_edges2()")
        dfile = util.read_dfile(filename, sep)
        result = []
        max_score = 0.0
        for line in dfile.lines:
            score = float(line[2])
            max_score = max(score, max_score)
            result.append((patches.patch_string_gene(organism_code, line[0]),
                           patches.patch_string_gene(organism_code, line[1]),
                           score))
        if not normalized:
            normalize_edges_to_max_score(result, max_score)
        return result

    def make_network(organism, ratios=None, check_size=False):
        """make network"""
        return network.Network.create("STRING", read_edges2(filename), weight,
                                      organism, ratios)

    return make_network


__all__ = ['get_network_factory2']
