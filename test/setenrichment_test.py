"""setenrichment_test.py - unit tests for set enrichment module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import json
import numpy as np
import cmonkey.set_enrichment as se
import cmonkey.datamatrix as dm
import cmonkey.util as util
import cmonkey.thesaurus as thesaurus

class DiscreteEnrichmentSetTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for DiscreteEnrichmentSet"""

    def test_construct(self):
        genes = {'gene1', 'gene2'}
        aset = se.DiscreteEnrichmentSet(genes)
        self.assertEquals(genes, aset.genes())
        self.assertEquals(genes, aset.genes_above_cutoff())


class CutoffEnrichmentSetTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for CutoffEnrichmentSet"""

    def test_construct(self):
        genes = {('gene1', 0.1), ('gene2', 0.6)}
        aset = se.CutoffEnrichmentSet(0.5, genes)
        self.assertEquals({'gene1', 'gene2'}, aset.genes())
        self.assertEquals({'gene2'}, aset.genes_above_cutoff())


class SetTypeTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for SetType"""

    def test_construct(self):
        genes1 = {'gene1', 'gene2'}
        aset1 = se.DiscreteEnrichmentSet(genes1)
        genes2 = {'gene3', 'gene4'}
        aset2 = se.DiscreteEnrichmentSet(genes2)
        set_type = se.SetType('mysettype', {'set1': aset1, 'set2': aset2}, 1.0)
        self.assertEquals('mysettype', set_type.name)
        self.assertEquals({'gene1', 'gene2', 'gene3', 'gene4'},
                          set_type.genes())


class MockMembership:
    def rows_for_cluster(self, cluster):
        return ['79073', '900', '29957', '11167', '4154', '84270', '84061', '4001', '26503', '4086', '51585', '4734', '5873']

class SetEnrichmentComputeClusterScoreTest(unittest.TestCase):

    def setUp(self):
        with open('testdata/hsa_mir_thesaurus.json') as infile:
            self.synonyms = json.load(infile)
        self.ratios = dm.create_from_csv('testdata/acc_rnaseq.tsv.gz', sep='\t')
        self.input_genes = self.ratios.row_names
        self.config_params = {'SetEnrichment': {'set_types': 'tfbs'},
                              'SetEnrichment-tfbs': {'set_file': 'testdata/test_sets.json', 'weight': 1.0} }
        self.set_types = se.read_set_types(self.config_params, self.synonyms, self.input_genes)
        self.canonical_rownames = set(map(lambda n: self.synonyms[n] if n in self.synonyms else n,
                                              self.ratios.row_names))

        self.canonical_row_indexes = {}
        for index, row in enumerate(self.ratios.row_names):
            if row in self.synonyms:
                self.canonical_row_indexes[self.synonyms[row]] = index
            else:
                self.canonical_row_indexes[row] = index

    def test_simple(self):
        num_clusters = 659
        cutoff = 0.05 / num_clusters
        ref_min_score = -4.702276
        scores, min_set, min_pvalue = se.compute_cluster_score_plain(1, cutoff, ref_min_score, self.ratios, MockMembership(),
                                                                     self.set_types[0], self.synonyms, self.canonical_rownames,
                                                                     self.canonical_row_indexes)
        self.assertEquals('hsa-miR-9', min_set)
        self.assertAlmostEquals(0.000212407251628, min_pvalue)
