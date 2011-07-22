"""cmonkey_test.py - unit test module for cmonkey module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from datamatrix import DataMatrix, DataMatrixCollection
from cmonkey import CMonkey, Membership

# setting up some simple dummy test input
# approximately 1 bicluster per 10 genes
NUM_GENES = 10
NUM_COLS = 2
GENE_NAMES = ['gene' + str(i) for i in range(NUM_GENES)]
COL_NAMES = ['cond' + str(i) for i in range(NUM_COLS)]
RATIOS = DataMatrix(NUM_GENES, NUM_COLS, GENE_NAMES, COL_NAMES)
RATIO_MATRICES = DataMatrixCollection([RATIOS])


class MockOrganism:
    """A mock organism"""

    def init_with(self, gene_names):
        pass

    def features(self):
        return {}

    def contigs(self):
        return []


class CMonkeyTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for CMonkey"""

    def test_create_cmonkey(self):
        """create a CMonkey object"""
        cmonkey = CMonkey(MockOrganism(), RATIO_MATRICES)
        self.assertFalse(cmonkey.finished())

    def test_create_cmonkey_with_config(self):
        """create CMonkey object, specifying an organism"""
        config = {'organism': 'homo sapiens'}
        cmonkey = CMonkey(MockOrganism(), RATIO_MATRICES, config)
        self.assertFalse(cmonkey.finished())

    def test_run_cmonkey_simple(self):
        """run CMonkey in the simplest way"""
        cmonkey = CMonkey(MockOrganism(), RATIO_MATRICES)
        cmonkey.run()
        self.assertTrue(cmonkey.finished())


class MembershipTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Membership"""

    def test_map_to_is_member_matrix(self):
        """tests the map_to_is_member_matrix function"""
        in_matrix = [[1, 2], [2, 3]]
        out = Membership.map_to_is_member_matrix(in_matrix, 3)
        self.assertEquals([[True, False], [True, True], [False, True]], out)


if __name__ == '__main__':
    unittest.main()
