"""meme_test.py - unit tests for meme module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import meme
import unittest


class MemeTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Location"""

    def test_read_meme_output(self):
        """tests the read_meme_output function"""
        with open('testdata/meme.out') as inputfile:
            motif_infos = meme.read_meme_output(inputfile)

        self.assertEquals(24, motif_infos[0].width())
        self.assertEquals(5, motif_infos[0].num_sites())
        self.assertEquals(95, motif_infos[0].llr())
        self.assertAlmostEquals(1.9e+3, motif_infos[0].evalue())
