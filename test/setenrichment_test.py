"""network_test.py - unit tests for network module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import set_enrichment as se

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
        set_type = se.SetType('mysettype', {'set1': aset1, 'set2': aset2})
        self.assertEquals('mysettype', set_type.name)
        self.assertEquals({'gene1', 'gene2', 'gene3', 'gene4'},
                          set_type.genes())
