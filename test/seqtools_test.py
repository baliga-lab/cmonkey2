"""seqtools_test.py - unit tests for seqtools module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from seqtools import subsequence, extract_upstream, subseq_counts
from seqtools import subseq_frequencies, markov_background
from seqtools import read_sequences_from_fasta_string


class SeqtoolsTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for seqtools."""

    def test_simple(self):
        """tests with a simple coordinate"""
        self.assertEquals('TTAG', subsequence("ATTAGCA", 2, 6))

    def test_reverse(self):
        """tests with simple coordinate, setting reverse flag"""
        self.assertEquals('CTAA', subsequence("ATTAGCA", 2, 6, reverse=True))

    def test_subseq_counts_1(self):
        """test subseq_counts() with length 1"""
        counts = subseq_counts(["ACCGTATA", "CACAT"], 1)
        self.assertEquals(5, counts['A'])
        self.assertEquals(4, counts['C'])
        self.assertEquals(1, counts['G'])
        self.assertEquals(3, counts['T'])

    def test_subseq_counts_2(self):
        """test subseq_counts() with length 2"""
        counts = subseq_counts(["ACCGTATA", "CACAT"], 2)
        self.assertEquals(2, counts['AC'])
        self.assertEquals(1, counts['CC'])
        self.assertEquals(1, counts['CG'])
        self.assertEquals(1, counts['GT'])
        self.assertEquals(2, counts['TA'])
        self.assertEquals(2, counts['AT'])
        self.assertEquals(2, counts['CA'])

    def test_subseq_frequencies_1(self):
        """test subseq_frequencies() with length 1"""
        freqs = subseq_frequencies(["ACCGTATA", "CACAT"], 1)
        self.assertAlmostEqual(5.0 / 13.0, freqs['A'])
        self.assertAlmostEqual(4.0 / 13.0, freqs['C'])
        self.assertAlmostEqual(1.0 / 13.0, freqs['G'])
        self.assertAlmostEqual(3.0 / 13.0, freqs['T'])

    def test_markov_background_0(self):
        """test markov_background() with order 0"""
        background = markov_background(["ACCGTATA", "CACAT"], 0)
        self.assertEquals(1, len(background))
        self.assertEquals(4, len(background[0]))

    def test_markov_background_1(self):
        """test markov_background() with order 1"""
        background = markov_background(["ACCGTATA", "CACAT"], 1)
        self.assertEquals(2, len(background))
        self.assertEquals(4, len(background[0]))
        self.assertEquals(7, len(background[1]))


class FastaTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for FASTA related functions"""

    def test_read_sequences_from_fasta_string(self):
        """test reading sequences from a string in FASTA format"""
        with open("testdata/fasta_test.fa") as inputfile:
            fasta_string = inputfile.read()
        seqs = read_sequences_from_fasta_string(fasta_string)
        print sorted(seqs.keys())
        self.assertEquals(7, len(seqs))
        seq = ("CCGAGGAAGACAGACGCAATTTCACATCGAACTCGTGTACGGCATCCTCT" +
               "TTATTGCCGGCTTTGCTTTTCTCGTCTTCCGCGTCGATCCCCGGGTGGCA" +
               "GCGTTCGAAGGAGGTCTCGTCATTGGTTACTTATTGAGAATTTAGGGGAA" +
               "AATGTCAATCTACGAGTGGA")
        self.assertEquals(seq, seqs['VNG6198H'])
