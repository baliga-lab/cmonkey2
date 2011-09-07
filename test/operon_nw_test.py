"""util_test.py - test classes for operon module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import seqtools as st
import organism as org
import microbes_online as mo
import util
import thesaurus as th


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
    def __init__(self, taxonomy_id, feature_map, synonyms=None):
        self.__taxonomy_id = taxonomy_id
        self.__feature_map = feature_map
        self.__synonyms = synonyms

    def taxonomy_id(self):
        return self.__taxonomy_id

    def features_for_genes(self, gene_names):
        return self.__feature_map

    def feature_id_for(self, gene_alias):
        return gene_alias


class MockOrganismWithSynonyms:
    """mock organism class with synonyms"""
    def __init__(self, taxonomy_id, feature_map, synonyms):
        self.__taxonomy_id = taxonomy_id
        self.__feature_map = feature_map
        self.__synonyms = synonyms

    def taxonomy_id(self):
        return self.__taxonomy_id

    def features_for_genes(self, gene_names):
        if self.__synonyms:
            result = {}
            for alias in gene_names:
                feature_id = self.__synonyms[alias]
                result[feature_id] = self.__feature_map[feature_id]
            return util.ThesaurusBasedMap(self.__synonyms, result)
        else:
            return self.__feature_map

    def feature_id_for(self, gene_alias):
        if self.__synonyms:
            return self.__synonyms[gene_alias]
        else:
            return gene_alias


class ReadOperonNetworkTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for read_from_microbes_online"""

    def test_build_operons(self):
        """tests building an operon list from two name lists"""
        names1 = ['VNG0001', 'VNG0007', 'VNG0008']
        names2 = ['VNG0003', 'VNG0008', 'VNG0009']
        operons = mo.build_operons(names1, names2)
        self.assertEquals([['VNG0001', 'VNG0003'],
                           ['VNG0007', 'VNG0008', 'VNG0009']], operons)

    def test_make_operon_edges_forward(self):
        """test when all genes of the operon are on the forward strand"""
        operon = ['gene1', 'gene2', 'gene3']
        features = {
            'gene1': st.Feature('feature1', 'typ1', 'feature_name1',
                                st.Location('contig1', 24, 89, False)),
            'gene2': st.Feature('feature2', 'typ1', 'feature_name2',
                                st.Location('contig1', 15, 21, False)),
            'gene3': st.Feature('feature3', 'typ2', 'feature_name3',
                                st.Location('contig1', 100, 154, False))
            }
        edges = mo.make_operon_pairs(operon, features)
        self.assertTrue(('gene2', 'gene1') in edges)
        self.assertTrue(('gene2', 'gene2') in edges)
        self.assertTrue(('gene2', 'gene3') in edges)
        self.assertEqual(3, len(edges))

    def test_make_operon_edges_reverse(self):
        """test when all genes of the operon are on the reverse strand"""
        operon = ['gene1', 'gene2', 'gene3']
        features = {
            'gene1': st.Feature('feature1', 'typ1', 'feature_name1',
                                st.Location('contig1', 24, 89, True)),
            'gene2': st.Feature('feature2', 'typ1', 'feature_name2',
                                st.Location('contig1', 15, 21, True)),
            'gene3': st.Feature('feature3', 'typ2', 'feature_name3',
                                st.Location('contig1', 100, 154, True))
            }
        edges = mo.make_operon_pairs(operon, features)
        self.assertTrue(('gene3', 'gene1') in edges)
        self.assertTrue(('gene3', 'gene2') in edges)
        self.assertTrue(('gene3', 'gene3') in edges)
        self.assertEqual(3, len(edges))

    def test_make_edges_from_predictions(self):
        """tests the make_edges_from_predictions function"""
        predictions = [('gene1', 'gene2'), ('gene2', 'gene3')]
        organism = MockOrganism('64091', {
                'gene1': st.Feature('feature1', 'typ1', 'feature_name1',
                                    st.Location('contig1', 24, 89, False)),
                'gene2': st.Feature('feature2', 'typ1', 'feature_name2',
                                    st.Location('contig1', 15, 21, False)),
                'gene3': st.Feature('feature3', 'typ2', 'feature_name3',
                                    st.Location('contig1', 100, 154, False))
                })
        edges = mo.make_pairs_from_predictions(predictions, organism)
        self.assertEquals([('gene2', 'gene1'), ('gene2', 'gene2'),
                           ('gene2', 'gene3')], edges)

    def test_get_network_factory(self):
        """test happy path"""
        microbes_online = MockMicrobesOnline('testdata/gnc64091.named')
        network = mo.get_network_factory(microbes_online, 20)(MockOrganism(
                '64091',
                 {'gene1': st.Feature('feature1', 'typ1', 'feature_name1',
                                      st.Location('contig1', 24, 89, False)),
                  'gene2': st.Feature('feature2', 'typ1', 'feature_name2',
                                      st.Location('contig1', 15, 21, False)),
                  'gene3': st.Feature('feature3', 'typ2', 'feature_name3',
                                      st.Location('contig1', 100, 154, False))
                  }))
        self.assertEquals(6, network.num_edges())
        self.assertEquals(6000, network.total_score())


class GetOperonPairsTest(unittest.TestCase):
    """integration test for operon pair retrieval"""

    def __make_organism(self):
        """makes a mock organism with almost real data"""
        features = {}
        dfile = util.DelimitedFile.read('testdata/Halobacterium_sp_features',
                                        comment='--')
        for line in dfile.lines():
            features[line[0]] = st.Feature(line[0], line[1], line[2],
                                           st.Location(line[3],
                                                       int(line[4]),
                                                       int(line[5]),
                                                       line[6] == 'R'))
        tfile = util.DelimitedFile.read(
            'testdata/Halobacterium_sp_feature_names', comment='--')
        synonyms = th.create_from_rsat_feature_names(tfile)
        return MockOrganismWithSynonyms('64091', features, synonyms)

    def __make_ref_operon_pairs(self):
        """returns reference operon pairs for comparison"""
        reffile = util.DelimitedFile.read(
            'testdata/operon_reftable.tsv', has_header=True, quote='"')
        refpairs = []
        for line in reffile.lines():
            refpairs.append((line[1], line[2]))
        return refpairs

    def test_make_operon_pairs(self):
        """test the make_operon_pairs() function in integration"""
        mo_db = MockMicrobesOnline('testdata/gnc64091_ref.named')
        organism = self.__make_organism()
        pairs = mo.get_operon_pairs(mo_db, organism)
        refpairs = self.__make_ref_operon_pairs()

        self.assertEquals(len(pairs), len(refpairs))
        for i in range(len(pairs)):
            self.assertEquals(refpairs[i], pairs[i])
