#!/usr/bin/env python3
"""meme_test.py - unit tests for meme module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import cmonkey.meme.meme as meme
import cmonkey.meme.mast as mast
import unittest


class MemeTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Location"""

    def test_read_meme430_output(self):
        """tests the read_meme_output function for MEME 4.3.0 format
        The MEME format essentially stayed the same from 4.3.0 until
        4.11.3. Starting from 4.11.4, motif naming convention changed.
        """
        with open('testdata/meme430.out') as inputfile:
            motif_infos = meme.from_text(inputfile.read(), 2)

        self.assertEqual(2, len(motif_infos))
        self.assertEqual(24, motif_infos[0].width)
        self.assertEqual(5, motif_infos[0].num_sites)
        self.assertEqual(95, motif_infos[0].llr)
        self.assertAlmostEqual(1.9e+3, motif_infos[0].evalue)
        self.assertEqual(5, len(motif_infos[0].sites))
        sites0 = motif_infos[0].sites
        self.assertEqual('NP_395673.1', sites0[0][0])
        self.assertEqual('NP_395728.1', sites0[1][0])
        self.assertEqual('NP_279261.1', sites0[2][0])
        self.assertEqual('NP_279811.1', sites0[3][0])
        self.assertEqual('NP_279211.1', sites0[4][0])

        self.assertEqual('+', sites0[0][1])
        self.assertEqual('-', sites0[1][1])
        self.assertEqual('-', sites0[2][1])
        self.assertEqual('+', sites0[3][1])
        self.assertEqual('-', sites0[4][1])

        self.assertEqual(8, sites0[0][2])
        self.assertEqual(138, sites0[1][2])
        self.assertEqual(27, sites0[2][2])
        self.assertEqual(4, sites0[3][2])
        self.assertEqual(7, sites0[4][2])

        self.assertAlmostEqual(2.61e-11, sites0[0][3])
        self.assertAlmostEqual(3.05e-10, sites0[1][3])
        self.assertAlmostEqual(9.59e-10, sites0[2][3])
        self.assertAlmostEqual(3.54e-09, sites0[3][3])
        self.assertAlmostEqual(5.81e-09, sites0[4][3])

        self.assertEqual('ACAGCGACAGCTTCCCGTCGATCT', sites0[0][5])
        self.assertEqual('AGATTGACATTTTCCCCTAAATTC', sites0[1][5])
        self.assertEqual('ACAGCAAAATCTACGTCTCGGACT', sites0[2][5])
        self.assertEqual('TGATAAAACACTTTATCTCTGTAT', sites0[3][5])
        self.assertEqual('ACGTAGACCGTATCGCGGAGATCT', sites0[4][5])

    def test_read_meme_4_11_4_output(self):
        """
        MEME 4.11.4 changes the MEME format in way that is non-compatible
        to how we handled the previous versions

        From the release notes:
          MEME motif names are now their consensus sequence, rather than a number,
          and the alternate ID is "MEME-i", where "i" = 1, 2, 3,... is the old motif
          name. MEME tests now check the integrity of all output formats: text, HTML and XML
        """
        with open('testdata/meme4.11.4.out') as inputfile:
            motif_infos = meme.from_text(inputfile.read(), 1)
            self.assertEqual(1, len(motif_infos))
            info = motif_infos[0]
            self.assertEqual(6, info.width)
            self.assertEqual(4, info.num_sites)
            self.assertEqual(35, info.llr)
            self.assertAlmostEqual(230.0, info.evalue)
            self.assertEqual(4, len(info.sites))

    def test_read_meme_4_12_0_output(self):
        """
        Just a test that the format reader has not changed in 4.12.0
        """
        with open('testdata/meme4.12.0.out') as inputfile:
            motif_infos = meme.from_text(inputfile.read(), 1)
            self.assertEqual(1, len(motif_infos))
            info = motif_infos[0]
            self.assertEqual(14, info.width)
            self.assertEqual(5, info.num_sites)
            self.assertEqual(68, info.llr)
            self.assertAlmostEqual(13.0, info.evalue)
            self.assertEqual(5, len(info.sites))

    def test_read_mast_output_430_1(self):
        """tests the read_mast_output_oldstyle function"""
        with open('testdata/mast430-1.out') as inputfile:
            pevalues, annotations = mast.from_430_text(
                inputfile.read(), ['VNG6198H', 'VNG0117H'])
        self.assertEqual('VNG6198H', pevalues[0][0])
        self.assertEqual('VNG0117H', pevalues[1][0])

        self.assertAlmostEqual(1.98e-14, pevalues[0][1])
        self.assertAlmostEqual(3.39e-12, pevalues[1][1])
        self.assertAlmostEqual(8.0e-12, pevalues[0][2])
        self.assertAlmostEqual(1.4e-09, pevalues[1][2])

        annot1 = list(annotations['VNG6198H'])
        self.assertAlmostEqual(6.6e-01, annot1[0][0])
        self.assertEqual(16, annot1[0][1])
        self.assertEqual(1, annot1[0][2])

        annot2 = list(annotations['VNG0117H'])
        self.assertAlmostEqual(2.0e-01, annot2[0][0])
        self.assertAlmostEqual(24, annot2[0][1])
        self.assertAlmostEqual(-1, annot2[0][2])
        self.assertAlmostEqual(4.9e-01, annot2[5][0])
        self.assertAlmostEqual(223, annot2[5][1])
        self.assertAlmostEqual(-2, annot2[5][2])

    def test_read_mast_output_430_2(self):
        """tests the read_mast_output function, this one has some
        more silly blank line placements"""
        with open('testdata/mast430-2.out') as inputfile:
            pevalues, annotations = mast.from_430_text(
                inputfile.read(), ['NP_279634.1', 'NP_279286.1'])
        self.assertTrue('NP_279634.1' in annotations)
        self.assertTrue('NP_279286.1' in annotations)

    def test_read_mast_output_430_3(self):
        """tests the read_mast_output function, this one has an incomplete block"""
        with open('testdata/mast430-3.out') as inputfile:
            pevalues, annotations = mast.from_430_text(
                inputfile.read(), ['NP_279608.1'])
        pev = [pevalue for pevalue in pevalues if pevalue[0] == 'NP_279608.1']
        self.assertAlmostEqual(3.9e-08, pev[0][1])
        self.assertAlmostEqual(9.61e-11, pev[0][2])
        self.assertTrue('NP_279608.1' in annotations)

    def test_read_mast_output_430_4(self):
        """tests the read_mast_output function, this has on sequence/annotation block"""
        with open('testdata/mast430-4.out') as inputfile:
            pevalues, annotations = mast.from_430_text(
                inputfile.read(), ['NP_280363.1', 'NP_280692.1'])
        pev = [pevalue for pevalue in pevalues if pevalue[0] == 'NP_280363.1']
        self.assertAlmostEqual(1.0, pev[0][1])
        self.assertAlmostEqual(4.0e02, pev[0][2])
        self.assertTrue('NP_280363.1' not in annotations)

    def test_read_mast_output_xml(self):
        """MEME 4.8.1 introduces XML based MAST files."""
        with open('testdata/mast-481.xml') as inputfile:
            pevalues, annotations = mast.from_xml_text(
                inputfile.read(), ['NP_280363.1', 'NP_280692.1'])
        pev = [pevalue for pevalue in pevalues if pevalue[0] == 'NP_280363.1']
        self.assertAlmostEqual(0.322, pev[0][1])
        self.assertAlmostEqual(130.0, pev[0][2])
        self.assertTrue('NP_280363.1' in annotations)

    def test_read_mast_output_4_11_xml(self):
        """MEME 4.11.x changed the XML format slightly"""
        with open('testdata/mast-4.11_output.xml') as inputfile:
            pevalues, annotations = mast.from_xml_text(
                inputfile.read(), ['NP_280363.1', 'NP_280692.1'])
        pev = [pevalue for pevalue in pevalues if pevalue[0] == 'NP_280363.1']
        self.assertAlmostEqual(1.0, pev[0][1])
        self.assertAlmostEqual(400.0, pev[0][2])
        self.assertTrue('NP_280363.1' in annotations)

    def test_read_mast_output_4_11_4_xml(self):
        """MAST 4.11.4 changed the naming according to MEME"""
        with open('testdata/mast-4.11.4_output.xml') as inputfile:
            pevalues, annotations = mast.from_xml_text(
                inputfile.read(), ['NP_280363.1', 'NP_280692.1'])
        pev = [pevalue for pevalue in pevalues if pevalue[0] == 'NP_280363.1']
        self.assertAlmostEqual(1.0, pev[0][1])
        self.assertAlmostEqual(400.0, pev[0][2])
        self.assertTrue('NP_280363.1' in annotations)


if __name__ == '__main__':
    unittest.main()
