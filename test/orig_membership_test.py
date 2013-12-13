"""orig_membership_test.py - unit test module for membership module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import membership as memb
import datamatrix as dm
import microarray as ma
import scoring

CONFIG_PARAMS = {
    'memb.clusters_per_row': 2,
    'memb.clusters_per_col': 5,
    'num_clusters': 43,
    'memb.prob_row_change': 0.5,
    'memb.prob_col_change': 1.0,
    'memb.max_changes_per_row': 1,
    'memb.max_changes_per_col': 5,
    'memb.min_cluster_rows_allowed': 1,
    'memb.max_cluster_rows_allowed': 5
}

class OrigMembershipTest(unittest.TestCase):
    """Test class for OrigMembership"""

    def test_constructor(self):
        """verify initialization"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3], 'C2': []},
                                CONFIG_PARAMS)
        self.assertEquals(43, m.num_clusters())
        self.assertEquals(2, m.num_clusters_per_row())
        self.assertEquals(5, m.num_clusters_per_column())
        self.assertEquals(0.5, m.probability_seeing_row_change())
        self.assertEquals(1.0, m.probability_seeing_col_change())
        self.assertEquals(1, m.max_changes_per_row())
        self.assertEquals(5, m.max_changes_per_col())
        self.assertEquals(1, m.min_cluster_rows_allowed())
        self.assertEquals(5, m.max_cluster_rows_allowed())
        self.assertEquals(0, m.min_cluster_columns_allowed())

        self.assertEquals([1, 5], m.row_memb['R1'].tolist())
        self.assertEquals([0, 0], m.row_memb['R2'].tolist())
        self.assertEquals([3, 0, 0, 0, 0], m.col_memb['C1'].tolist())
        self.assertEquals([0, 0, 0, 0, 0], m.col_memb['C2'].tolist())

        self.assertEquals({1, 5}, m.clusters_for_row('R1'))
        self.assertEquals({3}, m.clusters_for_column('C1'))

        self.assertEquals(2, m.num_clusters_for_row('R1'))
        self.assertEquals(1, m.num_clusters_for_column('C1'))

        self.assertEquals({'R1'}, m.rows_for_cluster(1))
        self.assertEquals({'C1'}, m.columns_for_cluster(3))
        self.assertEquals(1, m.num_row_members(1))
        self.assertEquals(0, m.num_row_members(7))
        self.assertEquals(1, m.num_column_members(3))
        self.assertEquals(0, m.num_column_members(7))

        self.assertEquals([2, 3], m.clusters_not_in_row('R1', [1, 2, 3, 5]))
        self.assertEquals([1, 2, 5], m.clusters_not_in_column('C1', [1, 2, 3, 5]))

        self.assertTrue(m.is_row_in_cluster('R1', 1))
        self.assertFalse(m.is_row_in_cluster('R1', 2))

        self.assertTrue(m.is_column_in_cluster('C1', 3))
        self.assertFalse(m.is_column_in_cluster('C1', 1))

    def test_add_cluster_to_row(self):
        """Happy path for add_cluster_to_row()"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3], 'C2': []},
                                CONFIG_PARAMS)
        m.add_cluster_to_row('R2', 3)
        self.assertEquals({3}, m.clusters_for_row('R2'))
        self.assertEquals(1, m.num_clusters_for_row('R2'))
        self.assertEquals({'R2'}, m.rows_for_cluster(3))

        m.add_cluster_to_row('R2', 1)
        self.assertEquals({1, 3}, m.clusters_for_row('R2'))
        self.assertEquals(2, m.num_clusters_for_row('R2'))
        self.assertEquals({'R1', 'R2'}, m.rows_for_cluster(1))
        self.assertEquals(2, m.num_row_members(1))

    def test_add_cluster_to_row_twice(self):
        """allow to add duplicate clusters"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3], 'C2': []},
                                CONFIG_PARAMS)
        m.add_cluster_to_row('R2', 1)
        m.add_cluster_to_row('R2', 1)
        self.assertEquals([1, 1], m.row_memb['R2'].tolist())
        self.assertEquals({1}, m.clusters_for_row('R2'))
        self.assertEquals(1, m.num_clusters_for_row('R2'))

    def test_add_cluster_to_row_error_full(self):
        """full row memberships can not be added to in the normal case"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3], 'C2': []},
                                CONFIG_PARAMS)
        self.assertRaises(Exception, m.add_cluster_to_row, 'R1', 2)

    def test_add_cluster_to_row_error_full_same_cluster(self):
        """add_cluster_to_row() will only work if there are enough 0-slots"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3], 'C2': []},
                                CONFIG_PARAMS)
        m.row_memb['R2'] = [2, 2]  # fill the spaces
        self.assertRaises(Exception, m.add_cluster_to_row, 'R2', 2)

    def test_add_cluster_to_row_full_force(self):
        """full row memberships can be added to with the force switch"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3], 'C2': []},
                                CONFIG_PARAMS)
        m.add_cluster_to_row('R1', 3, True)
        self.assertEquals({1, 3, 5}, m.clusters_for_row('R1'))
        self.assertEquals({'R1'}, m.rows_for_cluster(3))

    def test_add_cluster_to_column(self):
        """happy path for add_cluster_to_column()"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3], 'C2': []},
                                CONFIG_PARAMS)
        m.add_cluster_to_column('C1', 4)
        self.assertEquals({3, 4}, m.clusters_for_column('C1'))
        self.assertEquals(2, m.num_clusters_for_column('C1'))
        self.assertEquals({'C1'}, m.columns_for_cluster(4))

        m.add_cluster_to_column('C2', 3)
        self.assertEquals({'C1', 'C2'}, m.columns_for_cluster(3))
        self.assertEquals(2, m.num_column_members(3))

    def test_add_cluster_to_column_twice(self):
        """allow to add duplicate clusters"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3], 'C2': []},
                                CONFIG_PARAMS)
        m.add_cluster_to_column('C1', 4)
        m.add_cluster_to_column('C1', 4)
        self.assertEquals([3, 4, 4, 0, 0], m.col_memb['C1'].tolist())
        self.assertEquals({3, 4}, m.clusters_for_column('C1'))
        self.assertEquals(2, m.num_clusters_for_column('C1'))

    def test_add_cluster_to_column_error_full(self):
        """full row memberships can not be added to in the normal case"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3, 4, 5, 6, 7], 'C2': []},
                                CONFIG_PARAMS)
        self.assertRaises(Exception, m.add_cluster_to_column, 'C1', 2)

    def test_add_cluster_to_column_error_full_same_cluster(self):
        """full row memberships can not be added to in the normal case"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3], 'C2': []},
                                CONFIG_PARAMS)
        m.col_memb['C1'] = [1, 3, 6, 6, 6]
        self.assertRaises(Exception, m.add_cluster_to_column, 'C1', 2)

    def test_add_cluster_to_column_full_force(self):
        """full col memberships can be added to if forced"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3, 4, 5, 6, 7], 'C2': []},
                                CONFIG_PARAMS)
        m.add_cluster_to_column('C1', 2, force=True)
        self.assertEquals({2, 3, 4, 5, 6, 7}, m.clusters_for_column('C1'))
        self.assertEquals({'C1'}, m.columns_for_cluster(2))

    def test_replace_column_cluster(self):
        """full col memberships can be added to if forced"""
        m = memb.OrigMembership({'R1': [1, 5], 'R2': []}, {'C1': [3, 4, 3, 6, 7], 'C2': []},
                                CONFIG_PARAMS)
        m.replace_column_cluster('C1', 3, 1)
        self.assertEquals({1, 3, 4, 6, 7}, m.clusters_for_column('C1'))
        self.assertEquals({'C1'}, m.columns_for_cluster(1))


    def test_replace_row_cluster(self):
        """full col memberships can be added to if forced"""
        m = memb.OrigMembership({'R1': [1, 1], 'R2': []}, {'C1': [3, 4, 3, 6, 7], 'C2': []},
                                CONFIG_PARAMS)
        m.replace_row_cluster('R1', 1, 3)
        self.assertEquals({1, 3}, m.clusters_for_row('R1'))
        self.assertEquals({'R1'}, m.rows_for_cluster(3))

    def pickle_path(self):
        """returns the function-specific pickle-path"""
        return '%s/last_row_scores.pkl' % (self.__config_params['output_dir'])


if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(OrigMembershipTest))
    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
