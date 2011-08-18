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
            motif_infos = meme.read_meme_output(inputfile, 2)

        self.assertEquals(2, len(motif_infos))
        self.assertEquals(24, motif_infos[0].width())
        self.assertEquals(5, motif_infos[0].num_sites())
        self.assertEquals(95, motif_infos[0].llr())
        self.assertAlmostEquals(1.9e+3, motif_infos[0].evalue())
        self.assertEquals(5, len(motif_infos[0].sites()))
        sites0 = motif_infos[0].sites()
        self.assertEquals('NP_395673.1', sites0[0][0])
        self.assertEquals('NP_395728.1', sites0[1][0])
        self.assertEquals('NP_279261.1', sites0[2][0])
        self.assertEquals('NP_279811.1', sites0[3][0])
        self.assertEquals('NP_279211.1', sites0[4][0])

        self.assertEquals('+', sites0[0][1])
        self.assertEquals('-', sites0[1][1])
        self.assertEquals('-', sites0[2][1])
        self.assertEquals('+', sites0[3][1])
        self.assertEquals('-', sites0[4][1])

        self.assertEquals(8, sites0[0][2])
        self.assertEquals(138, sites0[1][2])
        self.assertEquals(27, sites0[2][2])
        self.assertEquals(4, sites0[3][2])
        self.assertEquals(7, sites0[4][2])

        self.assertAlmostEquals(2.61e-11, sites0[0][3])
        self.assertAlmostEquals(3.05e-10, sites0[1][3])
        self.assertAlmostEquals(9.59e-10, sites0[2][3])
        self.assertAlmostEquals(3.54e-09, sites0[3][3])
        self.assertAlmostEquals(5.81e-09, sites0[4][3])

        self.assertEquals('ACAGCGACAGCTTCCCGTCGATCT', sites0[0][4])
        self.assertEquals('AGATTGACATTTTCCCCTAAATTC', sites0[1][4])
        self.assertEquals('ACAGCAAAATCTACGTCTCGGACT', sites0[2][4])
        self.assertEquals('TGATAAAACACTTTATCTCTGTAT', sites0[3][4])
        self.assertEquals('ACGTAGACCGTATCGCGGAGATCT', sites0[4][4])
