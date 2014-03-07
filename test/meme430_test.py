"""meme430_test.py - integration tests for meme module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import meme
import motif
import unittest
import util
import rsat
import organism as org
import microbes_online
import stringdb
import datamatrix as dm
import os, os.path


KEGG_FILE_PATH = 'config/KEGG_taxonomy'
GO_FILE_PATH = 'config/proteome2taxid'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
CACHE_DIR = 'cache'

class FakeMembership:
    def num_clusters(self):
        return 1

    def rows_for_cluster(self, cluster):
        return ["VNG0457G", "VNG0715G", "VNG1143G", "VNG1182H", "VNG1190G",
                "VNG1408G", "VNG1551G", "VNG1562H", "VNG1698G", "VNG1760G",
                "VNG2191H", "VNG2199H", "VNG2344G", "VNG2410G", "VNG2567C"]


class Meme430Test(unittest.TestCase):  # pylint: disable-msg=R0904
    """This class tests a Halo setup"""

    def setUp(self):  # pylint: disable-msg=C0103
        if not os.path.exists('out'):
            os.mkdir('out')

    def test_meme_simple(self):
        """simplest of all: just run meme and parse the output, just tests
        if there will be appropriate output for the input"""
        meme_suite = meme.MemeSuite430()
        motif_infos, out = meme_suite.meme('testdata/meme_input1.fasta',
                                           'testdata/meme1.bg',
                                           num_motifs=1)
        self.assertEquals(1, len(motif_infos))
        self.assertEquals(24, motif_infos[0].width)
        self.assertEquals(3, motif_infos[0].num_sites)
        self.assertEquals(79, motif_infos[0].llr)
        self.assertAlmostEquals(1700, motif_infos[0].evalue)

    def test_motif_scoring(self):
        """tests the motif scoring in integration"""
        search_distances = {'upstream': (-20, 150)}
        scan_distances = {'upstream': (-30, 250)}

        matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
        infile = util.read_dfile('example_data/hal/halo_ratios5.tsv',
                                 has_header=True, quote='\"')
        ratio_matrix = matrix_factory.create_from(infile)
        organism = make_halo(ratio_matrix, search_distances, scan_distances)
        membership = FakeMembership()
        config_params = {'memb.min_cluster_rows_allowed': 3,
                         'memb.max_cluster_rows_allowed': 70,
                         'multiprocessing': False,
                         'num_clusters': 1,
                         'output_dir': 'out',
                         'debug': False,
                         'meme_version': '4.3.0',
                         'global_background': False,
                         'search_distances': {'upstream': (-20, 150)},
                         'nmotifs_rvec': 'c(rep(1, num_iterations/3), rep(2, num_iterations/3))',
                         'num_iterations': 2000,
                         'schedule': {'Motifs': lambda i: True},
                         'scaling': {'Motifs': ('scaling_const', 1.0)}}
        func = motif.MemeScoringFunction(organism, membership, ratio_matrix,
                                         motif_in_iteration=lambda x: True,
                                         config_params=config_params)
        iteration_result = { 'iteration': 100 }
        matrix = func.compute(iteration_result)

def make_halo(ratio_matrix, search_distances, scan_distances):
    """returns the organism object to work on"""
    keggfile = util.read_dfile(KEGG_FILE_PATH, comment='#')
    gofile = util.read_dfile(GO_FILE_PATH)
    rsatdb = rsat.RsatDatabase(rsat.RSAT_BASE_URL, CACHE_DIR)
    mo_db = microbes_online.MicrobesOnline(CACHE_DIR)
    stringfile = 'testdata/string_links_64091.tab'

    nw_factories = []
    if stringfile != None:
        nw_factories.append(stringdb.get_network_factory('hal', stringfile, 0.5))
    else:
        logging.warn("no STRING file specified !")

    nw_factories.append(microbes_online.get_network_factory(
            mo_db, max_operon_size=ratio_matrix.num_rows / 20, weight=0.5))

    org_factory = org.MicrobeFactory(org.make_kegg_code_mapper(keggfile),
                                     org.make_rsat_organism_mapper(rsatdb),
                                     org.make_go_taxonomy_mapper(gofile),
                                     mo_db,
                                     nw_factories)

    return org_factory.create('hal', search_distances, scan_distances)
