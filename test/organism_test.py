
"""organism_test.py - unit tests for organism module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import util
import network as nw
import organism as org
import seqtools as st


TAXONOMY_FILE_PATH = "testdata/KEGG_taxonomy"
PROT2TAXID_FILE_PATH = "testdata/proteome2taxid"
RSAT_LIST_FILE_PATH = "testdata/RSAT_genomes_listing.txt"

SEARCH_DISTANCES = {'upstream':(-20, 150)}
SCAN_DISTANCES = {'upstream':(-30, 250)}


# pylint: disable-msg=R0904
class KeggOrganismCodeMapperTest(unittest.TestCase):
    """Test class for KeggCodeMapper"""

    def test_get_existing_organism(self):
        """retrieve existing organism"""
        dfile = util.read_dfile(TAXONOMY_FILE_PATH, sep='\t',
                                has_header=True, comment='#')
        mapper = org.make_kegg_code_mapper(dfile)
        self.assertEquals('Helicobacter pylori 26695', mapper('hpy'))

    def test_get_non_existing_organism(self):
        """retrieve non-existing organism"""
        dfile = util.read_dfile(TAXONOMY_FILE_PATH, sep='\t',
                                has_header=True, comment='#')
        mapper = org.make_kegg_code_mapper(dfile)
        self.assertIsNone(mapper('nope'))


class GoTaxonomyMapperTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for get_go_taxonomy_id"""

    def test_get_existing(self):
        """retrieve an existing id"""
        dfile = util.read_dfile(PROT2TAXID_FILE_PATH, sep='\t',
                                has_header=False)
        mapper = org.make_go_taxonomy_mapper(dfile)
        self.assertEquals('64091', mapper('Halobacterium salinarium'))

    def test_get_non_existing(self):
        """retrieve None for a non-existing organism"""
        dfile = util.read_dfile(PROT2TAXID_FILE_PATH, sep='\t',
                                has_header=False)
        mapper = org.make_go_taxonomy_mapper(dfile)
        self.assertIsNone(mapper('does not exist'))


class MockRsatDatabase:
    """mock RsatDatabase"""

    def __init__(self, html):
        self.html = html

    def get_directory(self):
        """returns the directory listing's html text"""
        return self.html

    def get_organism(self, _):  # pylint: disable-msg=R0201
        """returns the organism file's content"""
        return '-- comment\nfoo bar; Eukaryota; something else'

    def get_organism_names(self, _):
        """returns a simulation of the organism_names.tab file"""
        return '-- comment\n4711\tRSAT organism'

    def get_features(self, _):
        """returns a fake feature.tab file"""
        return ('-- comment\n' +
                'NP_206803.1\tCDS\tnusB\tNC_000915.1\t123\t456\tD\n' +
                'NP_206804.1\tCDS\tnusC\tNC_000915.1\t234\t789\tR\n')

    def get_feature_names(self, _):
        """returns a fake feature_names.tab file"""
        return ("-- comment\nNP_206803.1\tNP_206803.1\tprimary\n" +
                "NP_206803.1\tVNG12345G\talternate\n" +
                "gene1\talt1\talternate")

    def get_contig_sequence(self, organism, contig):
        """return a contig sequence"""
        return "ACGTTTAAAAGAGAGAGAGACACAGTATATATTTTTTTAAAA"


class RsatOrganismMapperTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Tests the RsatOrganismMapper class"""

    def setUp(self):  # pylint: disable-msg=C0103
        """test fixture"""
        with open(RSAT_LIST_FILE_PATH) as inputfile:
            html = inputfile.read()
        self.mapper = org.make_rsat_organism_mapper(MockRsatDatabase(html))

    def test_mapper(self):
        """tests the get_organism method for an existing organism"""
        info = self.mapper('Halobacterium', None)
        self.assertEquals('Halobacterium_sp', info.species)
        self.assertTrue(info.is_eukaryote)
        self.assertEquals('4711', info.taxonomy_id)


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


def mock_go_mapper(rsat_organism):
    """A simple GO mock mapper to test whether the underscore is replaced
    in the factory"""
    if rsat_organism == 'RSAT organism':
        return 'GO taxonomy id'
    else:
        return None


class MockMicrobesOnline:
    def get_operon_predictions_for(self, organism_id):
        with open('testdata/gnc64091.named') as infile:
            return infile.read()


class MicrobeFactoryTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for OrganismFactory"""

    def test_create_prokaryote(self):
        """tests creating a Prokaryote"""
        factory = org.MicrobeFactory(
            lambda _: 'KEGG organism',
            lambda x, y, z: org.RsatSpeciesInfo(MockRsatDatabase(''),
                                                'RSAT_organism',
                                                False, 4711),
            mock_go_mapper,
            MockMicrobesOnline(),
            [])
        organism = factory.create('hpy', SEARCH_DISTANCES, SCAN_DISTANCES)
        self.assertEquals('hpy', organism.code)
        self.assertEquals('Hpy', organism.cog_organism())
        self.assertEquals('KEGG organism', organism.kegg_organism)
        self.assertEquals('RSAT_organism', organism.species())
        self.assertEquals('GO taxonomy id', organism.go_taxonomy_id)
        self.assertFalse(organism.is_eukaryote())
        self.assertIsNotNone(str(organism))

    def test_create_eukaryote(self):
        """tests creating an eukaryote"""
        factory = org.MicrobeFactory(
            lambda _: 'KEGG organism',
            lambda x, y, z: org.RsatSpeciesInfo(MockRsatDatabase(''),
                                             'RSAT_organism',
                                             True, 4711),
            lambda _: 'GO taxonomy id',
            MockMicrobesOnline(),
            [])
        organism = factory.create('hpy', SEARCH_DISTANCES, SCAN_DISTANCES)
        self.assertEquals('hpy', organism.code)
        self.assertTrue(organism.is_eukaryote())


class MicrobeTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Organism"""

    def setUp(self):  # pylint: disable-msg=C0103
        """test fixture"""
        self.mockFactory = MockNetworkFactory()
        self.organism = org.Microbe('hal', 'Halobacterium SP',
                                    org.RsatSpeciesInfo(MockRsatDatabase(''),
                                                        'Halobacterium_SP',
                                                        False,
                                                        12345),
                                    12345,
                                    MockMicrobesOnline(),
                                    [self.mockFactory],
                                    SEARCH_DISTANCES,
                                    SCAN_DISTANCES)

    def test_init_genome(self):
        """Tests the init_genome() method"""
        organism = self.organism
        scan_seqs = organism.sequences_for_genes_scan(['VNG12345G'], seqtype='upstream')
        self.assertEquals((st.Location('NC_000915.1', -128, 152, False),
                           'ACGTTTAAAAGAGAGAGAGACACAGTATATATTTTTTTAAAA'),
                          scan_seqs['NP_206803.1'])
        search_seqs = organism.sequences_for_genes_search(['VNG12345G'], seqtype='upstream')
        self.assertEquals((st.Location('NC_000915.1', -28, 142, False),
                           'ACGTTTAAAAGAGAGAGAGACACAGTATATATTTTTTTAAAA'),
                          search_seqs['NP_206803.1'])

    def test_search_no_operons(self):
        """Tests the init_genome() method"""
        organism = self.organism
        organism.use_operons = False
        scan_seqs = organism.sequences_for_genes_scan(['VNG12345G'], seqtype='upstream')
        self.assertEquals((st.Location('NC_000915.1', -128, 152, False),
                           'ACGTTTAAAAGAGAGAGAGACACAGTATATATTTTTTTAAAA'),
                          scan_seqs['NP_206803.1'])

    def test_get_networks(self):
        """tests the networks() method"""
        organism = self.organism
        networks = organism.networks()
        self.assertEquals(1, len(networks))
        self.assertEquals(organism, self.mockFactory.create_called_with)


class MockNetworkFactory:
    """a mock NetworkFactory"""
    def __init__(self):
        self.create_called_with = None

    def __call__(self, organism):
        self.create_called_with = organism
        return nw.Network.create('network', [], 0.0, check_size=False)
