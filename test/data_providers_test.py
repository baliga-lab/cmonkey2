"""Test classes for data_providers module"""
import unittest
from util import DelimitedFile
from data_providers import get_organism_for_code

TAXONOMY_FILE_PATH = "testdata/KEGG_taxonomy"


class OrganismCodeMappingTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for get_organism_for_code"""

    def test_get_existing_organism(self):
        """retrieve existing organism"""
        dfile = DelimitedFile.read(TAXONOMY_FILE_PATH, sep='\t',
                                   has_header=True, comment='#')
        self.assertEquals('Helicobacter pylori 26695',
                          get_organism_for_code(dfile, 'hpy'))

    def test_get_non_existing_organism(self):
        """retrieve non-existing organism"""
        dfile = DelimitedFile.read(TAXONOMY_FILE_PATH, sep='\t',
                                   has_header=True, comment='#')
        self.assertIsNone(get_organism_for_code(dfile, 'nope'))
