"""meme430_test.py - integration tests for meme module

This file is part of cm2meme. Please see README and LICENSE for
more information and licensing details.
"""
import cmonkey.meme_suite as meme
import unittest
import os, os.path


class Meme430Test(unittest.TestCase):  # pylint: disable-msg=R0904
    """This test class is a simple integration test with MEME 4.3.0
    It requires MEME 4.3.0 to be in your path and is not intended to
    run regularly.
    It simply tests if data is passed correctly.
    """
    def test_meme_simple(self):
        """simplest of all: just run meme and parse the output, just tests
        if there will be appropriate output for the input"""
        meme_suite = meme.MemeSuite430({'MEME': {'max_width': 24, 'background_order': 3,
                                                 'use_revcomp': 'True', 'arg_mod': 'zoops'}})
        motif_infos, out = meme_suite.meme('testdata/meme430_input1.fasta',
                                           'testdata/meme430_input1.bg',
                                           num_motifs=1)
        self.assertEqual(1, len(motif_infos))
        self.assertEqual(24, motif_infos[0].width)
        self.assertEqual(3, motif_infos[0].num_sites)
        self.assertEqual(79, motif_infos[0].llr)
        self.assertAlmostEqual(1700, motif_infos[0].evalue)


if __name__ == '__main__':
    unittest.main()
