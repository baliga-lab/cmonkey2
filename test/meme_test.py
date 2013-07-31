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
            motif_infos = meme.read_meme_output(inputfile.read(), 2)

        self.assertEquals(2, len(motif_infos))
        self.assertEquals(24, motif_infos[0].width)
        self.assertEquals(5, motif_infos[0].num_sites)
        self.assertEquals(95, motif_infos[0].llr)
        self.assertAlmostEquals(1.9e+3, motif_infos[0].evalue)
        self.assertEquals(5, len(motif_infos[0].sites))
        sites0 = motif_infos[0].sites
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

    def test_read_mast_output_oldstyle(self):
        """tests the read_mast_output function"""
        with open('testdata/mast.out') as inputfile:
            pevalues, annotations = meme.read_mast_output_oldstyle(
                inputfile.read(), ['VNG6198H', 'VNG0117H'])
        self.assertEquals('VNG6198H', pevalues[0][0])
        self.assertEquals('VNG0117H', pevalues[1][0])

        self.assertAlmostEquals(1.98e-14, pevalues[0][1])
        self.assertAlmostEquals(3.39e-12, pevalues[1][1])
        self.assertAlmostEquals(8.0e-12, pevalues[0][2])
        self.assertAlmostEquals(1.4e-09, pevalues[1][2])

        annot1 = annotations['VNG6198H']
        self.assertAlmostEquals(6.6e-01, annot1[0][0])
        self.assertEquals(16, annot1[0][1])
        self.assertEquals(1, annot1[0][2])

        annot2 = annotations['VNG0117H']
        self.assertAlmostEquals(2.0e-01, annot2[0][0])
        self.assertAlmostEquals(24, annot2[0][1])
        self.assertAlmostEquals(-1, annot2[0][2])
        self.assertAlmostEquals(4.9e-01, annot2[5][0])
        self.assertAlmostEquals(223, annot2[5][1])
        self.assertAlmostEquals(-2, annot2[5][2])

    def test_read_mast_output_oldstyle2(self):
        """tests the read_mast_output function, this one has some
        more silly blank line placements"""
        with open('testdata/mast2.out') as inputfile:
            pevalues, annotations = meme.read_mast_output_oldstyle(
                inputfile.read(), ['NP_279634.1', 'NP_279286.1'])
        self.assertTrue('NP_279634.1' in annotations)
        self.assertTrue('NP_279286.1' in annotations)

    def test_read_mast_output_oldstyle3(self):
        """tests the read_mast_output function, this one has an incomplete block"""
        with open('testdata/mast3.out') as inputfile:
            pevalues, annotations = meme.read_mast_output_oldstyle(
                inputfile.read(), ['NP_279608.1'])
        pev = [pevalue for pevalue in pevalues if pevalue[0] == 'NP_279608.1']
        self.assertAlmostEquals(3.9e-08, pev[0][1])
        self.assertAlmostEquals(9.61e-11, pev[0][2])
        self.assertTrue('NP_279608.1' in annotations)

    def test_read_mast_output_oldstyle4(self):
        """tests the read_mast_output function, this has on sequence/annotation block"""
        with open('testdata/mast4.out') as inputfile:
            pevalues, annotations = meme.read_mast_output_oldstyle(
                inputfile.read(), ['NP_280363.1', 'NP_280692.1'])
        pev = [pevalue for pevalue in pevalues if pevalue[0] == 'NP_280363.1']
        self.assertAlmostEquals(1.0, pev[0][1])
        self.assertAlmostEquals(4.0e02, pev[0][2])
        self.assertTrue('NP_280363.1' not in annotations)

    def test_read_mast_output_xml(self):
        with open('testdata/mast-481.xml') as inputfile:
            pevalues, annotations = meme.read_mast_output_xml(
                inputfile.read(), ['NP_280363.1', 'NP_280692.1'])
        pev = [pevalue for pevalue in pevalues if pevalue[0] == 'NP_280363.1']
        self.assertAlmostEquals(0.322, pev[0][1])
        self.assertAlmostEquals(130.0, pev[0][2])
        self.assertTrue('NP_280363.1' in annotations)
