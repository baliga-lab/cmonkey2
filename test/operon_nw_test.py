"""util_test.py - test classes for operon module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from microbes_online import get_network_factory, build_operons
from microbes_online import make_operon_edges, make_edges_from_predictions
from organism import Feature


class MockMicrobesOnline:
    """mock class for MicrobesOnline"""

    def __init__(self, filename):
        """creates a mock instance"""
        with open(filename) as inputfile:
            self.content = inputfile.read()

    def get_operon_predictions_for(self, _):
        """mocked MicrobesOnline operon prediction result"""
        return self.content


class MockOrganism:
    """mock Organism class"""
    def __init__(self, taxonomy_id, feature_map):
        self.__taxonomy_id = taxonomy_id
        self.__feature_map = feature_map

    def taxonomy_id(self):
        return self.__taxonomy_id

    def get_feature(self, gene):
        return self.__feature_map[gene]


class ReadOperonNetworkTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for read_from_microbes_online"""

    def test_build_operons(self):
        """tests building an operon list from two name lists"""
        names1 = ['VNG0001', 'VNG0007', 'VNG0008']
        names2 = ['VNG0003', 'VNG0008', 'VNG0009']
        operons = build_operons(names1, names2)
        self.assertEquals([['VNG0001', 'VNG0003'],
                           ['VNG0007', 'VNG0008', 'VNG0009']], operons)

    def test_make_operon_edges_forward(self):
        """test when all genes of the operon are on the forward strand"""
        operon = ['gene1', 'gene2', 'gene3']
        organism = MockOrganism('64091', {
                'gene1': Feature('feature1', 'typ1', 'feature_name1',
                                  'contig1', 24, 89, False),
                'gene2': Feature('feature2', 'typ1', 'feature_name2',
                                  'contig1', 15, 21, False),
                'gene3': Feature('feature3', 'typ2', 'feature_name3',
                                  'contig1', 100, 154, False)
                })
        edges = make_operon_edges(operon, organism)
        self.assertTrue(('gene2', 'gene1') in edges)
        self.assertTrue(('gene2', 'gene2') in edges)
        self.assertTrue(('gene2', 'gene3') in edges)
        self.assertEqual(3, len(edges))

    def test_make_operon_edges_reverse(self):
        """test when all genes of the operon are on the reverse strand"""
        operon = ['gene1', 'gene2', 'gene3']
        organism = MockOrganism('64091', {
                'gene1': Feature('feature1', 'typ1', 'feature_name1',
                                  'contig1', 24, 89, True),
                'gene2': Feature('feature2', 'typ1', 'feature_name2',
                                  'contig1', 15, 21, True),
                'gene3': Feature('feature3', 'typ2', 'feature_name3',
                                  'contig1', 100, 154, True)
                })
        edges = make_operon_edges(operon, organism)
        self.assertTrue(('gene3', 'gene1') in edges)
        self.assertTrue(('gene3', 'gene2') in edges)
        self.assertTrue(('gene3', 'gene3') in edges)
        self.assertEqual(3, len(edges))

    def test_make_edges_from_predictions(self):
        """tests the make_edges_from_predictions function"""
        predictions = [('gene1', 'gene2'), ('gene2', 'gene3')]
        organism = MockOrganism('64091', {
                'gene1': Feature('feature1', 'typ1', 'feature_name1',
                                  'contig1', 24, 89, False),
                'gene2': Feature('feature2', 'typ1', 'feature_name2',
                                  'contig1', 15, 21, False),
                'gene3': Feature('feature3', 'typ2', 'feature_name3',
                                  'contig1', 100, 154, False)
                })
        edges = make_edges_from_predictions(predictions, organism)
        self.assertEquals([('gene2', 'gene1'), ('gene2', 'gene2'),
                           ('gene2', 'gene3')], edges)

    def test_get_network_factory(self):
        """test happy path"""
        microbes_online = MockMicrobesOnline('testdata/gnc64091.named')
        network = get_network_factory(microbes_online)(MockOrganism(
                '64091',
                 {'gene1': Feature('feature1', 'typ1', 'feature_name1',
                                    'contig1', 24, 89, False),
                  'gene2': Feature('feature2', 'typ1', 'feature_name2',
                                    'contig1', 15, 21, False),
                  'gene3': Feature('feature3', 'typ2', 'feature_name3',
                                    'contig1', 100, 154, False)}))
        self.assertEquals(5, network.num_edges())
        self.assertEquals(5000, network.total_score())
