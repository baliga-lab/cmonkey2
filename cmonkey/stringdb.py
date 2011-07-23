"""stringdb.py - cMonkey STRING database interface module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
from util import DelimitedFile
import logging
import re
import math
from network import NetworkEdge, Network

STRING_FILE2 = 'string_links_64091.tab'
PROTEIN_PREFIX = re.compile('^string:\d+[.]')


def read_edges(filename):
    """read the edges from a tab separated file. This is the original
    STRING format"""
    logging.info("stringdb.read_edges()")
    dfile = DelimitedFile.read(filename)
    result = []
    exists = {}
    max_score = 0.0
    for line in dfile.lines():
        protein1 = re.sub(PROTEIN_PREFIX, '', line[0])
        protein2 = re.sub(PROTEIN_PREFIX, '', line[1])
        combined_score = line[14].split('|')[0].split(':')[1]

        key = ':'.join([protein1, protein2, combined_score])
        if not key in exists:
            exists[key] = True
            score = float(combined_score)
            result.append(NetworkEdge(protein1, protein2, score))
            if score > max_score:
                max_score = score

    # normalize scores to 1000
    for edge in result:
        score = edge.score() / max_score * 1000.0
        score = 1000 * math.exp(score / 1000.0) / math.exp(1.0)
        edge.set_score(score)
    return result


def read_edges2(filename):
    """just read a preprocessed file, much faster to debug"""
    logging.info("stringdb.read_edges2()")
    dfile = DelimitedFile.read(filename)
    result = []
    for line in dfile.lines():
        result.append(NetworkEdge(line[0], line[1], float(line[2])))
    return result


def get_network_factory(filename):
    """temporary STRING network factory"""
    def make_network(_):
        """make network"""
        return Network.create("STRING", read_edges2(filename))

    return make_network

__all__ = ['get_network_factory']
