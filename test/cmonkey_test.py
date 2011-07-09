"""cmonkey_test.py - unit test module for cmonkey module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from datatypes import DataMatrix, DataMatrixCollection
from cmonkey import CMonkey, Membership

# setting up some simple dummy test input
# approximately 1 bicluster per 10 genes
NUM_GENES = 10
NUM_COLS = 2
GENE_NAMES = ['gene' + str(i) for i in range(NUM_GENES)]
COL_NAMES = ['cond' + str(i) for i in range(NUM_COLS)]
RATIOS = DataMatrix(NUM_GENES, NUM_COLS, GENE_NAMES, COL_NAMES)
RATIO_MATRICES = DataMatrixCollection([RATIOS])
DUMMY_ORGANISM = None


class CMonkeyTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for CMonkey"""

    def test_create_cmonkey(self):
        """create a CMonkey object"""
        cmonkey = CMonkey(DUMMY_ORGANISM, RATIO_MATRICES)
        self.assertFalse(cmonkey.run_finished)
        self.assertEquals('hpy', cmonkey.configuration['organism'])

    def test_create_cmonkey_with_config(self):
        """create CMonkey object, specifying an organism"""
        config = {'organism': 'homo sapiens'}
        cmonkey = CMonkey(DUMMY_ORGANISM, RATIO_MATRICES, config)
        self.assertFalse(cmonkey.run_finished)
        self.assertEquals('homo sapiens', cmonkey.configuration['organism'])

    def test_run_cmonkey_simple(self):
        """run CMonkey in the simplest way"""
        cmonkey = CMonkey(DUMMY_ORGANISM, RATIO_MATRICES)
        cmonkey.run()
        self.assertTrue(cmonkey.run_finished)


class MembershipTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Membership"""

    def test_map_to_is_member_matrix(self):
        """tests the map_to_is_member_matrix function"""
        in_matrix = [[1, 2], [2, 3]]
        out = Membership.map_to_is_member_matrix(in_matrix, 3)
        self.assertEquals([[True, False], [True, True], [False, True]], out)


if __name__ == '__main__':
    unittest.main()
