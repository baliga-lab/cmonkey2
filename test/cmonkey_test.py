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
    def __init__(self, networks):
        self.__networks = networks

    def init_with(self, gene_names):
        pass

    def features(self):
        return {}

    def contigs(self):
        return []

    def networks(self):
        return self.__networks

class MockNetwork:
    """A mock network"""
    def __init__(self, total_score):
        self.__total = total_score

    def total_score(self):
        return self.__total

    def normalize_scores_to(self, score):
        self.__total = score


class CMonkeyTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for CMonkey"""

    def test_create_cmonkey(self):
        """create a CMonkey object"""
        cmonkey = CMonkey(MockOrganism([]), RATIO_MATRICES)
        self.assertFalse(cmonkey.finished())

    def test_run_cmonkey_simple(self):
        """run CMonkey in the simplest way"""
        cmonkey = CMonkey(MockOrganism([]), RATIO_MATRICES)
        cmonkey.run()
        self.assertTrue(cmonkey.finished())

    def test_run_cmonkey_with_networks(self):
        """run CMonkey in the simplest way"""
        network1 = MockNetwork(3)
        network2 = MockNetwork(5)
        cmonkey = CMonkey(MockOrganism([network1, network2]), RATIO_MATRICES)
        cmonkey.run()
        self.assertTrue(cmonkey.finished())
        self.assertEquals(5, network1.total_score())
        self.assertEquals(5, network2.total_score())


class MembershipTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Membership"""

    def test_map_to_is_member_matrix(self):
        """tests the map_to_is_member_matrix function"""
        in_matrix = [[1, 2], [2, 3]]
        out = Membership.map_to_is_member_matrix(in_matrix, 3)
        self.assertEquals([[True, False], [True, True], [False, True]], out)


if __name__ == '__main__':
    unittest.main()
