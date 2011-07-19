"""seqtools_test.py - unit tests for seqtools module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from seqtools import subsequence, extract_upstream


class SeqtoolsTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for seqtools."""

    def test_simple(self):
        """tests with a simple coordinate"""
        self.assertEquals('TTAG', subsequence("ATTAGCA", 2, 6))

    def test_reverse(self):
        """tests with simple coordinate, setting reverse flag"""
        self.assertEquals('CTAA', subsequence("ATTAGCA", 2, 6, reverse=True))
