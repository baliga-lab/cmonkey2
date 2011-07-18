"""thesaurus_test.py - test classes for thesaurus module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from thesaurus import DelimitedFileThesaurusFactory


class MockDelimitedFile:
    """just a mocked DelimitedFile"""

    def lines(self):
        """mocked lines() method"""
        return [['alt1', 'gene1'], ['alt2', 'gene1'], ['alt3', 'gene2']]


class DelimitedFileFactoryTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for DelimitedFileThesaurusFactory"""

    def test_create(self):
        factory = DelimitedFileThesaurusFactory(MockDelimitedFile())
        thesaurus = factory.create()
        self.assertIsNotNone(thesaurus)
        self.assertEquals('gene1', thesaurus['alt1'])
        self.assertEquals('gene1', thesaurus['alt2'])
        self.assertEquals('gene2', thesaurus['alt3'])
