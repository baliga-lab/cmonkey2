"""seqtools_test.py - unit tests for seqtools module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import os
import seqtools as st


class LocationTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Location"""

    def test_equals_true(self):
        """Tests the equality operator"""
        location1 = st.Location("contig", 1, 2, True)
        location2 = st.Location("contig", 1, 2, True)

        self.assertTrue(location1 == location1)
        self.assertTrue(location1 == location2)

    def test_equals_false(self):
        """Tests the equality operator"""
        location1 = st.Location("contig1", 1, 2, True)
        location2 = st.Location("contig1", 1, 2, False)
        location3 = st.Location("contig2", 1, 2, True)

        self.assertFalse(location1 == location2)
        self.assertFalse(location1 == location3)

    def test_notequals(self):
        """Tests the equality operator"""
        location1 = st.Location("contig1", 1, 2, True)
        location2 = st.Location("contig1", 1, 2, False)
        location3 = st.Location("contig2", 1, 2, True)

        self.assertTrue(location1 != location2)
        self.assertTrue(location1 != location3)
        self.assertFalse(location1 != location1)


class SeqtoolsTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for seqtools."""

    def test_simple(self):
        """tests with a simple coordinate"""
        self.assertEquals('TTAG', st.subsequence("ATTAGCA", 2, 6))

    def test_reverse(self):
        """tests with simple coordinate, setting reverse flag"""
        self.assertEquals('CTAA', st.subsequence("ATTAGCA", 2, 6, reverse=True))

    def test_subseq_counts_1(self):
        """test subseq_counts() with length 1"""
        counts = st.subseq_counts(["ACCGTATA", "CACAT"], 1)
        self.assertEquals(5, counts['A'])
        self.assertEquals(4, counts['C'])
        self.assertEquals(1, counts['G'])
        self.assertEquals(3, counts['T'])

    def test_subseq_counts_2(self):
        """test subseq_counts() with length 2"""
        counts = st.subseq_counts(["ACCGTATA", "CACAT"], 2)
        self.assertEquals(2, counts['AC'])
        self.assertEquals(1, counts['CC'])
        self.assertEquals(1, counts['CG'])
        self.assertEquals(1, counts['GT'])
        self.assertEquals(2, counts['TA'])
        self.assertEquals(2, counts['AT'])
        self.assertEquals(2, counts['CA'])

    def test_subseq_frequencies_1(self):
        """test subseq_frequencies() with length 1"""
        freqs = st.subseq_frequencies(["ACCGTATA", "CACAT"], 1)
        self.assertAlmostEqual(5.0 / 13.0, freqs['A'])
        self.assertAlmostEqual(4.0 / 13.0, freqs['C'])
        self.assertAlmostEqual(1.0 / 13.0, freqs['G'])
        self.assertAlmostEqual(3.0 / 13.0, freqs['T'])

    def test_markov_background_0(self):
        """test markov_background() with order 0"""
        background = st.markov_background(["ACCGTATA", "CACAT"], 0)
        self.assertEquals(1, len(background))
        self.assertEquals(4, len(background[0]))

    def test_markov_background_1(self):
        """test markov_background() with order 1"""
        background = st.markov_background(["ACCGTATA", "CACAT"], 1)
        self.assertEquals(2, len(background))
        self.assertEquals(4, len(background[0]))
        self.assertEquals(7, len(background[1]))


class FastaTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for FASTA related functions"""

    def tearDown(self):
        """cleanup"""
        if os.path.exists('/tmp/fasta_tmp.fa'):
            os.remove('/tmp/fasta_tmp.fa')

    def test_read_sequences_from_fasta_string(self):
        """test reading sequences from a string in FASTA format"""
        with open("testdata/fasta_test.fa") as inputfile:
            fasta_string = inputfile.read()
        seqs = st.read_sequences_from_fasta_string(fasta_string)
        self.assertEquals(7, len(seqs))
        seq = ("CCGAGGAAGACAGACGCAATTTCACATCGAACTCGTGTACGGCATCCTCT" +
               "TTATTGCCGGCTTTGCTTTTCTCGTCTTCCGCGTCGATCCCCGGGTGGCA" +
               "GCGTTCGAAGGAGGTCTCGTCATTGGTTACTTATTGAGAATTTAGGGGAA" +
               "AATGTCAATCTACGAGTGGA")
        self.assertEquals('VNG6198H', seqs[6][0])
        self.assertEquals(seq, seqs[6][1])

    def test_read_sequences_from_fasta_file(self):
        """test reading sequences from a string in FASTA format"""
        with open("testdata/fasta_test.fa") as inputfile:
            fasta_string = inputfile.read()
        seqs = st.read_sequences_from_fasta_file('testdata/fasta_test.fa')
        self.assertEquals(7, len(seqs))
        seq = ("CCGAGGAAGACAGACGCAATTTCACATCGAACTCGTGTACGGCATCCTCT" +
               "TTATTGCCGGCTTTGCTTTTCTCGTCTTCCGCGTCGATCCCCGGGTGGCA" +
               "GCGTTCGAAGGAGGTCTCGTCATTGGTTACTTATTGAGAATTTAGGGGAA" +
               "AATGTCAATCTACGAGTGGA")
        self.assertEquals('VNG6198H', seqs[6][0])
        self.assertEquals(seq, seqs[6][1])

    def test_write_sequences_to_fasta_file(self):
        """Tests writing to a FASTA file"""
        seqs = st.read_sequences_from_fasta_file('testdata/fasta_test.fa')
        st.write_sequences_to_fasta_file(seqs, '/tmp/fasta_tmp.fa')
        seqs2 = st.read_sequences_from_fasta_file('/tmp/fasta_tmp.fa')
        self.assertEquals(seqs, seqs2)
