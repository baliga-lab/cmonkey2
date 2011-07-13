"""organism_test.py - unit tests for organism module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from util import DelimitedFile
from organism import make_kegg_code_mapper, make_go_taxonomy_mapper
from organism import make_rsat_organism_mapper, OrganismFactory

TAXONOMY_FILE_PATH = "testdata/KEGG_taxonomy"
PROT2TAXID_FILE_PATH = "testdata/proteome2taxid"
RSAT_LIST_FILE_PATH = "testdata/RSAT_genomes_listing.txt"


# pylint: disable-msg=R0904
class KeggOrganismCodeMapperTest(unittest.TestCase):
    """Test class for KeggCodeMapper"""

    def test_get_existing_organism(self):
        """retrieve existing organism"""
        dfile = DelimitedFile.read(TAXONOMY_FILE_PATH, sep='\t',
                                   has_header=True, comment='#')
        mapper = make_kegg_code_mapper(dfile)
        self.assertEquals('Helicobacter pylori 26695', mapper('hpy'))

    def test_get_non_existing_organism(self):
        """retrieve non-existing organism"""
        dfile = DelimitedFile.read(TAXONOMY_FILE_PATH, sep='\t',
                                   has_header=True, comment='#')
        mapper = make_kegg_code_mapper(dfile)
        self.assertIsNone(mapper('nope'))


class GoTaxonomyMapperTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for get_go_taxonomy_id"""

    def test_get_existing(self):
        """retrieve an existing id"""
        dfile = DelimitedFile.read(PROT2TAXID_FILE_PATH, sep='\t',
                                   has_header=False)
        mapper = make_go_taxonomy_mapper(dfile)
        self.assertEquals('64091', mapper('Halobacterium salinarium'))

    def test_get_non_existing(self):
        """retrieve None for a non-existing organism"""
        dfile = DelimitedFile.read(PROT2TAXID_FILE_PATH, sep='\t',
                                   has_header=False)
        mapper = make_go_taxonomy_mapper(dfile)
        self.assertIsNone(mapper('does not exist'))


class MockRsatDatabase:
    """mock RsatDatabase"""

    def __init__(self, html):
        self.html = html

    def get_directory_html(self):
        """returns the directory listing's html text"""
        return self.html

    def get_organism_file(self, _):  # pylint: disable-msg=R0201
        """returns the organism file's content"""
        return 'foo bar; Eukaryota; something else'


class RsatOrganismMapperTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Tests the RsatOrganismMapper class"""

    def setUp(self):  # pylint: disable-msg=C0103
        """test fixture"""
        with open(RSAT_LIST_FILE_PATH) as inputfile:
            html = inputfile.read()
        self.mapper = make_rsat_organism_mapper(MockRsatDatabase(html))

    def test_mapper(self):
        """tests the get_organism method for an existing organism"""
        self.assertEquals(('Halobacterium_sp', True),
                          self.mapper('Halobacterium'))


class MockRsatOrganismMapper:
    """mock RSAT organism mapper"""

    def __init__(self, is_eukaryotic):  # pylint: disable-msg=R0201
        """create an instance of this mock mapper"""
        self.is_eukaryotic = is_eukaryotic

    def get_organism(self, _):  # pylint: disable-msg=R0201
        """returns an organism for a KEGG organism"""
        return "RSAT organism"

    def is_eukaryote(self, _):
        """determine whether eukaryotic or prokaryotic"""
        return self.is_eukaryotic


class OrganismTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Organism"""

    def test_create_prokaryote(self):
        """tests creating a Prokaryote"""
        factory = OrganismFactory(lambda _: 'KEGG organism',
                                  lambda _: ('RSAT organism', False),
                                  lambda _: 'GO taxonomy id')
        organism = factory.create('hpy')
        self.assertEquals('hpy', organism.code)
        self.assertEquals('Hpy', organism.get_cog_organism())
        self.assertEquals('KEGG organism', organism.kegg_organism)
        self.assertEquals('RSAT organism', organism.rsat_organism)
        self.assertEquals('GO taxonomy id', organism.go_taxonomy_id)
        self.assertFalse(organism.is_eukaryote())

    def test_create_eukaryote(self):
        """tests creating an eukaryote"""
        factory = OrganismFactory(lambda _: 'KEGG organism',
                                  lambda _: ('RSAT organism', True),
                                  lambda _: 'GO taxonomy id')
        organism = factory.create('hpy')
        self.assertEquals('hpy', organism.code)
        self.assertTrue(organism.is_eukaryote())
