"""thesaurus_test.py - test classes for thesaurus module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from thesaurus import DelimitedFileThesaurusFactory1
from thesaurus import DelimitedFileThesaurusFactory2


class MockDelimitedFile1:
    """just a mocked DelimitedFile"""

    def lines(self):
        """mocked lines() method"""
        return [['alt1', 'gene1'], ['alt2', 'gene1'], ['alt3', 'gene2']]


class MockDelimitedFile2:
    """just a mocked DelimitedFile"""

    def lines(self):
        """mocked lines() method"""
        return [['gene1', 'alt1;alt2'], ['gene2', 'alt3']]


class DelimitedFileFactoryTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for DelimitedFileThesaurusFactory"""

    def test_factory1(self):
        """test the delimited file factory simple version"""
        factory = DelimitedFileThesaurusFactory1(MockDelimitedFile1())
        thesaurus = factory.create()
        self.assertIsNotNone(thesaurus)
        self.assertEquals('gene1', thesaurus['alt1'])
        self.assertEquals('gene1', thesaurus['alt2'])
        self.assertEquals('gene2', thesaurus['alt3'])

    def test_factory1(self):
        """test the delimited file factory second version"""
        factory = DelimitedFileThesaurusFactory2(MockDelimitedFile2())
        thesaurus = factory.create()
        self.assertIsNotNone(thesaurus)
        self.assertEquals('gene1', thesaurus['alt1'])
        self.assertEquals('gene1', thesaurus['alt2'])
        self.assertEquals('gene2', thesaurus['alt3'])
