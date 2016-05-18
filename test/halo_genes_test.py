"""halo_genes_test.py - integration tests for halo organism

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import cmonkey.organism as org
import cmonkey.microbes_online as microbes_online
import unittest
import cmonkey.util as util
import cmonkey.rsat as rsat
import testutil

class HaloGeneTest(unittest.TestCase):
    """Tests for retrieving gene sequences in Halo"""
    def test_get_one_sequence(self):
        """get one simple sequence"""
        search_distances = {'upstream': (-20, 150)}
        scan_distances = {'upstream': (-30, 250)}
        halo = testutil.make_halo(search_distances, scan_distances)
        print("------")
        print("VNG1551G: ", halo.features_for_genes(['VNG1551G']))
        print("VNG1550G: ", halo.features_for_genes(['VNG1550G']))
        print("VNG1561C: ", halo.features_for_genes(['VNG1561C']))
        print("------")
        seq = halo.sequences_for_genes_search(['VNG1551G'], 'upstream')
        self.assertTrue('NP_280354.1' in seq)
        print("FINAL LOC: ", seq['NP_280354.1'][0])
        """
        operon_map = halo.operon_map()
        for key, value in operon_map.items():
            print("'%s' -> '%s'" % (key, value))
        print("OPERON MAPPED TO: ", halo.operon_map()['NP_280354.1'])
        """
        self.assertEquals('GTGATTCGACCATTACTGCAAGTTCAGACGACCCCAATTCAAGTAGTTTTGTGTAACCGCCGGCGTCGGGGGCGCTCGCGCCCATCTAAGAAAGCTCACTTTCCCTAATACAATCAAAATTGTTTTGGGTGCTTCTGACGTTGTGCCACCGATGGCACAGACACAGCTCC',
                          seq['NP_280354.1'][1])
