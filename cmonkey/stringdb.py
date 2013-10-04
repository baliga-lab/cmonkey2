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


def read_edges(filename, sep='\t', prefix=PROTEIN_PREFIX):
    """read the edges from a tab separated file. This is the original
    STRING format used in R cMonkey
    """
    logging.info("stringdb.read_edges()")
    dfile = util.DelimitedFile.read(filename, sep)
    result = []
    exists = {}
    max_score = 0.0
    for line in dfile.lines():
        protein1 = re.sub(prefix, '', line[0])
        protein2 = re.sub(prefix, '', line[1])
        combined_score = line[14].split('|')[0].split(':')[1]

        key = ':'.join([protein1, protein2, combined_score])
        if not key in exists:
            exists[key] = True
            score = float(combined_score)
            result.append([protein1, protein2, score])
            if score > max_score:
                max_score = score
    return normalize_edges_to_max_score(result, max_score)


def normalize_edges_to_max_score(edges, max_score):
    """normalize scores to 1000, for combined scores"""
    for edge in edges:
        score = edge[2] / max_score * 1000.0
        score = 1000 * math.exp(score / 1000.0) / math.exp(1.0)
        edge[2] = score
    return edges


def get_network_factory(filename, weight):
    """temporary STRING network factory using the default STRING format
    This is the legacy factory method from R/cMonkey which supports
    the old STRING database format
    """
    def make_network(_):
        """make network"""
        return network.Network.create("STRING",
                                      read_edges(filename), weight)

    return make_network


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
            result.append([patches.patch_string_gene(organism_code, line[0]),
                           patches.patch_string_gene(organism_code, line[1]),
                           score])
        if not normalized:
            normalize_edges_to_max_score(result, max_score)
        return result

    def make_network(_):
        """make network"""
        return network.Network.create("STRING", read_edges2(filename), weight)

    return make_network


def get_network_factory2_FS(filename, weight, sep='\t'):
    """STRING network factory from preprocessed edge file
    (protein1, protein2, combined_score), scores are already
    normalized to 1000"""
    def read_edges2(filename):
        """just read a preprocessed file, much faster to debug"""
        logging.info("\x1b[31mstringdb:\t\x1b[0mreading interaction network - stringdb.read_edges2()")
        dfile = util.read_dfile(filename, sep)
        result = []
        for line in dfile.lines:
            result.append((line[0], line[1], float(line[2])))
        return result

    def make_network(_):
        """make network"""
        return network.FS_Network("STRING", read_edges2(filename), weight)

    return make_network


def get_network_factory3(filename, weight):
    """STRING network factory from preprocessed edge file
    (row, protein1, protein2, combined_score), scores are not yet
    normalized to 1000.
    This is a customized factory method used for certain human setups.
    """
    def read_edges3(filename):
        """just read a preprocessed file, much faster to debug"""
        logging.info("stringdb.read_edges3()")
        dfile = util.read_dfile(filename, sep=",", has_header=True, quote='"')
        result = []
        for line in dfile.lines:
            result.append([line[1], line[2], float(line[3])])
        return result

    def make_network(_):
        """make network"""
        return network.Network.create("STRING", read_edges3(filename), weight)

    return make_network


__all__ = ['get_network_factory', 'get_network_factory2',
           'get_network_factory3']
