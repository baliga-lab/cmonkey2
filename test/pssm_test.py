"""pssm_test.py - test classes for pssm module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import cmonkey.pssm as p

class PssmTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Pssm"""

    def setUp(self):  # pylint: disable-msg=C0103
        """test fixture"""
        pass

    def test_create_empty(self):
        """Test creation with no data"""
        pssm = p.Pssm('pssm')
        self.assertEquals('pssm', pssm.name)
        self.assertIsNone(pssm.sites)
        self.assertIsNone(pssm.e_value)

    def test_read_fasta(self):
        """tests reading Pssm's from a FASTA file"""
        pssms = []
        with open('testdata/three_miRNAs.fasta') as infile:
            pssms = p.read_fasta(infile)
        self.assertEquals(3, len(pssms))
        self.assertEquals('hsa-miR-1', pssms[0].name)
        self.assertEquals('hsa-miR-7', pssms[1].name)
        self.assertEquals('hsa-miR-124', pssms[2].name)

        self.assertAlmostEquals(0.80, pssms[0][0][0])
        self.assertAlmostEquals(0.04, pssms[0][0][1])
        self.assertAlmostEquals(0.05, pssms[0][0][2])
        self.assertAlmostEquals(0.11, pssms[0][0][3])

        self.assertAlmostEquals(0.02, pssms[1][1][0])
        self.assertAlmostEquals(0.01, pssms[1][1][1])
        self.assertAlmostEquals(0.02, pssms[1][1][2])
        self.assertAlmostEquals(0.85, pssms[1][1][3])

        self.assertAlmostEquals(0.10, pssms[2][2][0])
        self.assertAlmostEquals(0.10, pssms[2][2][1])
        self.assertAlmostEquals(0.75, pssms[2][2][2])
        self.assertAlmostEquals(0.05, pssms[2][2][3])

        self.assertEquals('ACATTCCA', pssms[0].consensus_motif())
        self.assertEquals('GTCTTCCA', pssms[1].consensus_motif())
        self.assertEquals('GTGCCTNA', pssms[2].consensus_motif())

        self.assertEquals('ACATTCCA', pssms[0].consensus_motif(three=False))
        self.assertEquals('GTCTTCCA', pssms[1].consensus_motif(three=False))
        self.assertEquals('GTGCCTNA', pssms[2].consensus_motif(three=False))
