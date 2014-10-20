"""meme430_test.py - integration tests for meme module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import meme
import motif
import unittest
import util
import organism as org
import datamatrix as dm
import os, os.path

import testutil

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
        meme_suite = meme.MemeSuite430({'MEME': {'max_width': 24, 'background_order': 3,
                                                 'use_revcomp': 'True', 'arg_mod': 'zoops'}})
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
        organism = testutil.make_halo(search_distances, scan_distances, ratio_matrix)
        membership = FakeMembership()
        config_params = {'memb.min_cluster_rows_allowed': 3,
                         'memb.max_cluster_rows_allowed': 70,
                         'multiprocessing': False,
                         'num_clusters': 1,
                         'output_dir': 'out',
                         'debug': {},
                         'search_distances': {'upstream': (-20, 150)},
                         'num_iterations': 2000,
                         'MEME': {'schedule': lambda i: True,
                                  'version': '4.3.0',
                                  'global_background': False,
                                  'arg_mod': 'zoops',
                                  'nmotifs_rvec': 'c(rep(1, num_iterations/3), rep(2, num_iterations/3))',
                                  'use_revcomp': 'True', 'max_width': 24, 'background_order': 3},
                         'Motifs': {'schedule': lambda i: True, 'scaling': ('scaling_const', 1.0)}}
        func = motif.MemeScoringFunction(organism, membership, ratio_matrix,
                                         config_params=config_params)
        iteration_result = { 'iteration': 100 }
        matrix = func.compute(iteration_result)

