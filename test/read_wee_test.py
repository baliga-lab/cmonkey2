"""read_wee_test.py - test classes for wee module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import weeder as w

class ReadWeeTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Pssm"""

    def test_read(self):
        """Test creation with no data"""
        reader = w.WeederReader('testdata/perm_miR1_0_1.fasta.wee', 'base')
        reader.read()
        self.assertEquals(('GGGGGC', 1.7), reader.top_hit6())
        self.assertEquals(('CTGGGGGG', 1.89), reader.top_hit8())
        seqnames = reader.sequence_names()
        for num in range(len(seqnames)):
            self.assertEquals('perm' + str(num), seqnames[num])
        self.assertEquals(2, len(reader.pssms()))

        self.assertEquals('base_GGGGGC', reader.pssms()[0].name)
        self.assertEquals(1.7, reader.pssms()[0].e_value)
        self.assertEquals(137, len(reader.pssms()[0].sites))

        self.assertEquals('base_CTGGGGGG', reader.pssms()[1].name)
        self.assertEquals(1.89, reader.pssms()[1].e_value)
        self.assertEquals(25, len(reader.pssms()[1].sites))
