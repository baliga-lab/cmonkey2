"""pssm_test.py - test classes for pssm module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import datamatrix as dm
import scoring as s

class CombinerTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Pssm"""

    def test_combine_single(self):
        """Test combine with a single matrix"""
        m = dm.DataMatrix(2, 2, [[0.1, 0.2], [0.1, 0.2]])
        result = s.combine([m], [1.0], None, True)

