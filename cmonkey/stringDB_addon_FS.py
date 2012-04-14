'''
Created on Apr 2, 2012

@author: frank
'''

import logging
import re
import math
import util
import network_addon_FS as nw_FS


def get_network_factory2_FS(filename, weight, sep='\t'):
    """STRING network factory from preprocessed edge file
    (protein1, protein2, combined_score), scores are already
    normalized to 1000"""
    def read_edges2(filename):
        """just read a preprocessed file, much faster to debug"""
        logging.info("\x1b[31mstringdb:\t\x1b[0mreading interaction network - stringdb.read_edges2()")
        dfile = util.DelimitedFile.read(filename, sep)
        result = []
        for line in dfile.lines():
            result.append((line[0], line[1],
                                              float(line[2])))
        logging.info("\x1b[31mstringdb:\t\x1b[0mreading interaction network - done")
        return result

    def make_network(_):
        """make network"""
        tnw = nw_FS.Network("STRING", weight)
        tnw.buildNWfromList(read_edges2(filename))
        return tnw
    
    return make_network

def Babel_network_factory1_FS(NWclass, weight):
    """STRING network factory from preprocessed edge file
    (protein1, protein2, combined_score), scores are already
    normalized to 1000"""
    def getNW(NWclass):
        """just read a preprocessed file, much faster to debug"""
        logging.info("\x1b[31mstringdb:\t\x1b[0mreading interaction network - BabelFish")
        NWclass.EstablishCommunication()
        NWclass.ImportSTRING()
        NWclass.KillCommunication()
        logging.info("\x1b[31mstringdb:\t\x1b[0mreading interaction network - done")
        return NWclass.getNetwork()

    def make_network(_):
        """make network"""
        tnw = nw_FS.Network("STRING", weight)
        tnw.putGraph(getNW(NWclass))
        return tnw
    
    return make_network



