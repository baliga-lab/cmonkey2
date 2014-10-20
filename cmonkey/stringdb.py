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


def get_network_factory(organism_code, filename, weight, sep='\t',
                        normalized=False):
    """STRING network factory from preprocessed edge file
    (protein1, protein2, combined_score), scores are already
    normalized to 1000.
    This is the standard factory method used for Microbes.
    """
    def can_add_edge(node1, node2, thesaurus, cano_genes):
        """check whether we can add the edge"""
        if cano_genes is not None:
            return (node1 in thesaurus and node2 in thesaurus
                    and thesaurus[node1] in cano_genes and thesaurus[node2] in cano_genes)
        else:
            return node1 in thesaurus and node2 in thesaurus

    def read_edges2(filename, organism, ratios):
        """just read a preprocessed file, much faster to debug"""
        logging.info("stringdb.read_edges2()")
        dfile = util.read_dfile(filename, sep)
        result = []
        max_score = 0.0
        thesaurus = organism.thesaurus()
        if ratios:
            cano_genes = {thesaurus[row] for row in ratios.row_names
                          if row in thesaurus}
        else:
            cano_genes = None

        num_ignored = 0

        for line in dfile.lines:
            node1 = patches.patch_string_gene(organism_code, line[0])
            node2 = patches.patch_string_gene(organism_code, line[1])
            score = float(line[2])
            max_score = max(score, max_score)

            if can_add_edge(node1, node2, thesaurus, cano_genes):
                result.append((intern(node1), intern(node2), score))
            else:
                num_ignored += 1

        if not normalized:
            result = normalize_edges_to_max_score(result, max_score)

        logging.info("stringdb.read_edges2(), %d edges read, %d edges ignored",
                     len(result), num_ignored)
        return result

    def make_network(organism, ratios=None, check_size=False):
        """make network"""
        return network.Network.create("STRING",
                                      read_edges2(filename, organism, ratios),
                                      weight,
                                      organism, ratios)

    return make_network


__all__ = ['get_network_factory']
