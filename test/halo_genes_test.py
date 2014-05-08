"""halo_genes_test.py - integration tests for halo organism

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import organism as org
import microbes_online
import unittest
import util
import rsat

KEGG_FILE_PATH = 'config/KEGG_taxonomy'
GO_FILE_PATH = 'config/proteome2taxid'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
CACHE_DIR = 'cache'

class HaloGeneTest(unittest.TestCase):
    """Tests for retrieving gene sequences in Halo"""
    def test_get_one_sequence(self):
        """get one simple sequence"""
        search_distances = {'upstream': (-20, 150)}
        scan_distances = {'upstream': (-30, 250)}
        halo = make_halo(search_distances, scan_distances)
        print "------"
        print "VNG1551G: ", halo.features_for_genes(['VNG1551G'])
        print "VNG1550G: ", halo.features_for_genes(['VNG1550G'])
        print "VNG1561C: ", halo.features_for_genes(['VNG1561C'])
        print "------"
        seq = halo.sequences_for_genes_search(['VNG1551G'], 'upstream')
        self.assertTrue('NP_280354.1' in seq)
        print "FINAL LOC: ", seq['NP_280354.1'][0]
        """
        operon_map = halo.operon_map()
        for key, value in operon_map.items():
            print "'%s' -> '%s'" % (key, value)
        print "OPERON MAPPED TO: ", halo.operon_map()['NP_280354.1']
        """
        self.assertEquals('GTGATTCGACCATTACTGCAAGTTCAGACGACCCCAATTCAAGTAGTTTTGTGTAACCGCCGGCGTCGGGGGCGCTCGCGCCCATCTAAGAAAGCTCACTTTCCCTAATACAATCAAAATTGTTTTGGGTGCTTCTGACGTTGTGCCACCGATGGCACAGACACAGCTCC',
                          seq['NP_280354.1'][1])

    
def make_halo(search_distances, scan_distances):
    """returns the organism object to work on"""
    keggfile = util.read_dfile(KEGG_FILE_PATH, comment='#')
    gofile = util.read_dfile(GO_FILE_PATH)
    rsatdb = rsat.RsatDatabase(rsat.RSAT_BASE_URL, CACHE_DIR)
    mo_db = microbes_online.MicrobesOnline(CACHE_DIR)

    keggorg = org.make_kegg_code_mapper(keggfile)('hal')
    rsat_info = org.RsatSpeciesInfo(rsatdb, keggorg, None, None)
    gotax = org.make_go_taxonomy_mapper(gofile)(rsat_info.go_species())
    return org.Microbe('hal', keggorg, rsat_info, gotax, mo_db, [],
                       search_distances, scan_distances, True, None)
