"""thesaurus_test.py - test classes for thesaurus module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import thesaurus


class MockDelimitedFile1:
    """just a mocked DelimitedFile"""

    def __init__(self):
        """mocked lines() method"""
        self.lines = [['alt1', 'gene1'], ['alt2', 'gene1'], ['alt3', 'gene2']]


class MockDelimitedFile2:
    """just a mocked DelimitedFile"""

    def __init__(self):
        """mocked lines() method"""
        self.lines = [['gene1', 'alt1;alt2'], ['gene2', 'alt3']]


class MockRsatFeatureNameFile:
    """a mocked RSAT feature names delimited file"""

    def __init__(self):
        """mocked lines() method"""
        self.lines = [['NAME1', 'PRIME1', 'primary'], ['NAME1', 'ALT1', 'alternate'],
                      ['NAME2', 'PRIME2', 'primary'], ['NAME2', 'VNG2664Gm']]


class DelimitedFileFactoryTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for DelimitedFileThesaurusFactory"""

    def test_create_from_delimited_file1(self):
        """test the delimited file simple version"""
        thes = thesaurus.create_from_delimited_file1(MockDelimitedFile1())
        self.assertEquals('gene1', thes['alt1'])
        self.assertEquals('gene1', thes['alt2'])
        self.assertEquals('gene2', thes['alt3'])

    def test_create_from_delimited_file2(self):
        """test the delimited file second version"""
        thes = thesaurus.create_from_delimited_file2(MockDelimitedFile2())
        self.assertEquals('GENE1', thes['ALT1'])
        self.assertEquals('GENE1', thes['ALT2'])
        self.assertEquals('GENE2', thes['ALT3'])

    def test_create_from_rsat_feature_names_no_transform(self):
        """test the creation from RSAT feature names file"""
        thes = thesaurus.create_from_rsat_feature_names(
            MockRsatFeatureNameFile())
        self.assertEquals('NAME1', thes['PRIME1'])
        self.assertEquals('NAME1', thes['ALT1'])
        self.assertEquals('NAME2', thes['PRIME2'])
        self.assertEquals('NAME2', thes['VNG2664Gm'])

    def test_create_from_rsat_feature_names_with_transform(self):
        """test the creation from RSAT feature names using a key transformer"""
        thes = thesaurus.create_from_rsat_feature_names(
            MockRsatFeatureNameFile(), [lambda x: [x, x.rstrip('m')]])
        self.assertEquals('NAME1', thes['PRIME1'])
        self.assertEquals('NAME1', thes['ALT1'])
        self.assertEquals('NAME2', thes['PRIME2'])
        self.assertEquals('NAME2', thes['VNG2664G'])
