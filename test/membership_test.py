"""membership_test.py - unit test module for membership module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import membership as memb
import datamatrix as dm
import microarray as ma
import scoring


class MockSeedRowMemberships:

    def __init__(self):
        self.was_called = False

    def __call__(self, row_membership, data_matrix):
        self.was_called = True


class MockSeedColumnMemberships:

    def __init__(self):
        self.was_called = False

    def __call__(self, data_matrix, row_membership, num_clusters,
                 num_clusters_per_column):
        self.was_called = True
        return [[0], [0], [0], [0], [0]]

MATRIX1 = [[-0.8682966, 0.0, -0.34731863, 1.786210, 0.0],
           [-0.7642219, 0.008938267, -0.33071589, 1.850221, 0.0],
           [-0.5507068, 0.045892230, -0.22028270, 1.991723, 0.0],
           [-0.3798518, 0.061836346, -0.24734538, 2.058267, 0.0],
           [-0.5656293, 0.087864752, -0.09335630, 2.020889, 0.0],
           [-0.6259898, 0.071135199, -0.14227040, 1.977559, 0.0],
           [-0.6856569, 0.0, -0.23374668, 1.924514, 0.07791556],
           [-0.5205586, 0.0, -0.09943255, 2.041292, 0.08773460],
           [-0.8318400, 0.175124204, -0.06254436, 1.882585, 0.0],
           [-0.6615701, 0.0, -0.06343823, 1.975648, 0.07250084]]


class MockDataMatrix:

    def __init__(self, num_rows):
        self.num_rows = num_rows
        self.num_columns = 5
        self.row_names = [('GENE%d' % index) for index in range(self.num_rows)]
        self.column_names = [('COND%d' % index) for index in range(self.num_columns)]
        self.values = MATRIX1

    def submatrix(self, rows, cols):
        return self


class ClusterMembershipTest(unittest.TestCase):
    """Test class for ClusterMembership"""

    def test_constructor(self):
        """tests the constructor"""
        membership = memb.ClusterMembership(
            {'R1': [1, 3], 'R2': [2, 3]},
            {'C1': [1, 2], 'C2': [2]},
            {'num_clusters': 3,
             'memb.clusters_per_row': 2})
        self.assertEquals([1, 3], membership.clusters_for_row('R1'))
        self.assertEquals([2, 3], membership.clusters_for_row('R2'))
        self.assertEquals([], membership.clusters_for_row('R3'))
        self.assertEquals([1, 2], membership.clusters_for_column('C1'))
        self.assertEquals([2], membership.clusters_for_column('C2'))
        self.assertEquals([], membership.clusters_for_column('C3'))
        self.assertEquals(['R1'], membership.rows_for_cluster(1))
        self.assertEquals(['R2'], membership.rows_for_cluster(2))
        self.assertEquals(['R1', 'R2'], membership.rows_for_cluster(3))
        self.assertEquals(['C1'], membership.columns_for_cluster(1))
        self.assertEquals(3, membership.num_clusters())
        self.assertEquals(2, membership.num_clusters_per_row())
        self.assertEquals(1, membership.num_row_members(1))
        self.assertEquals(2, membership.num_row_members(3))
        self.assertEquals(0, membership.num_row_members(4))
        self.assertEquals(1, membership.num_column_members(1))
        self.assertEquals(2, membership.num_column_members(2))
        self.assertEquals(0, membership.num_column_members(3))

    def test_create_membership(self):
        """test creating a ClusterMembership object"""
        datamatrix = MockDataMatrix(3)
        seed_row_memberships = MockSeedRowMemberships()
        seed_col_memberships = MockSeedColumnMemberships()
        membership = memb.ClusterMembership.create(datamatrix,
                                                   seed_row_memberships,
                                                   seed_col_memberships,
                                                   {'memb.clusters_per_row': 2,
                                                    'num_clusters': 1,
                                                    'memb.clusters_per_col': 2})
        self.assertTrue(seed_row_memberships.was_called)
        self.assertTrue(seed_col_memberships.was_called)
        self.assertEquals(1, membership.num_clusters())

    def test_create_membership_use_defaults(self):
        """test creating a ClusterMembership object using default parameters"""
        datamatrix = MockDataMatrix(3)
        seed_row_memberships = MockSeedRowMemberships()
        seed_col_memberships = MockSeedColumnMemberships()
        membership = memb.ClusterMembership.create(datamatrix,
                                                   seed_row_memberships,
                                                   seed_col_memberships,
                                                   {'memb.clusters_per_col': 2,
                                                    'memb.clusters_per_row': 2,
                                                    'num_clusters': 43})
        self.assertTrue(seed_row_memberships.was_called)
        self.assertTrue(seed_col_memberships.was_called)
        self.assertEquals(memb.NUM_CLUSTERS, membership.num_clusters())

    def test_compute_column_scores_submatrix(self):
        """tests compute_column_scores_submatrix"""
        matrix = dm.DataMatrix(10, 5, ['R1', 'R2', 'R3', 'R4', 'R5', 'R6',
                                       'R7', 'R8', 'R9', 'R10'],
                               ['C1', 'C2', 'C3', 'C4', 'C5'],
                               MATRIX1)
        result = scoring.compute_column_scores_submatrix(matrix)
        scores = result.values[0]
        self.assertEqual(5, len(scores))
        self.assertAlmostEqual(0.03085775, scores[0])
        self.assertAlmostEqual(0.05290099, scores[1])
        self.assertAlmostEqual(0.05277032, scores[2])
        self.assertAlmostEqual(0.00358045, scores[3])
        self.assertAlmostEqual(0.03948821, scores[4])

    def test_seed_column_members(self):
        """tests seed_column_members"""
        data_matrix = dm.DataMatrix(3, 2, ["GENE1", "GENE2", "GENE3"],
                                    ['COL1', 'COL2'], [[1.0, 2.0], [2.0, 1.0],
                                                    [2.0, 1.0]])
        row_membership = [[1, 0], [2, 0], [1, 0]]
        column_members = ma.seed_column_members(data_matrix, row_membership,
                                                2, 2)
        self.assertEquals([2, 1], column_members[0])
        self.assertEquals([2, 1], column_members[1])

    def test_add_cluster_to_row(self):
        """tests adding a cluster to a row"""
        membership = memb.ClusterMembership(
            row_is_member_of={'R1': [1, 3], 'R2': [2, 3]},
            column_is_member_of={},
            config_params={'memb.clusters_per_col': 5,
                           'memb.clusters_per_row': 2})
        membership.add_cluster_to_row('R3', 1)
        self.assertEquals([1], membership.clusters_for_row('R3'))
        self.assertEquals(['R1', 'R3'], membership.rows_for_cluster(1))

    def test_add_cluster_to_row_exceed_limit(self):
        """tests adding a row to a cluster, checking the limit"""
        membership = memb.ClusterMembership(
            row_is_member_of={'R1': [1, 3], 'R2': [2, 3]},
            column_is_member_of={},
            config_params={})
        self.assertRaises(Exception,
                          membership.add_cluster_to_row, 'R1', 2)

    def test_add_cluster_to_column(self):
        """tests adding a cluster to a column"""
        membership = memb.ClusterMembership(
            row_is_member_of={},
            column_is_member_of={'C1': [1, 3], 'C2': [2, 3]},
            config_params={'memb.clusters_per_col': 2,
                           'memb.clusters_per_row': 1})
        membership.add_cluster_to_column('C3', 1)
        membership.add_cluster_to_column('C3', 1)
        self.assertEquals([1], membership.clusters_for_column('C3'))
        self.assertEquals(['C1', 'C3'], membership.columns_for_cluster(1))

    def test_add_cluster_to_column_twice(self):
        """tests adding a cluster to a column twice and making sure it's only in there once"""
        membership = memb.ClusterMembership(
            row_is_member_of={},
            column_is_member_of={'C1': [1, 3], 'C2': [2, 3]},
            config_params={'memb.clusters_per_col': 2,
                           'memb.clusters_per_row': 1})
        membership.add_cluster_to_column('C3', 1)
        self.assertEquals([1], membership.clusters_for_column('C3'))
        self.assertEquals(['C1', 'C3'], membership.columns_for_cluster(1))


    def test_add_cluster_to_column_exceed_limit(self):
        """tests adding a column to a cluster, hitting the limit"""
        membership = memb.ClusterMembership(
            row_is_member_of={},
            column_is_member_of={'C1': [1, 3], 'C2': [2, 3]},
            config_params={})
        self.assertRaises(Exception,
                          membership.add_cluster_to_column, 'C1', 2)

    def test_remove_cluster_from_row(self):
        """tests removing a row from a cluster"""
        membership = memb.ClusterMembership(
            row_is_member_of={'R1': [1, 3], 'R2': [2, 3]},
            column_is_member_of={},
            config_params={})
        membership.remove_cluster_from_row('R1', 3)
        self.assertEquals([1], membership.clusters_for_row('R1'))
        self.assertEquals(['R2'], membership.rows_for_cluster(3))

    def test_remove_cluster_from_column(self):
        """tests removing a column from a cluster"""
        membership = memb.ClusterMembership(
            row_is_member_of={},
            column_is_member_of={'C1': [1, 3], 'C2': [2, 3]},
            config_params={})
        membership.remove_cluster_from_column('C1', 3)
        self.assertEquals([1], membership.clusters_for_column('C1'))
        self.assertEquals(['C2'], membership.columns_for_cluster(3))

    def test_replace_column_cluster(self):
        membership = memb.ClusterMembership(
            row_is_member_of={},
            column_is_member_of={'C1': [1, 3], 'C2': [2, 3]},
            config_params={'memb.clusters_per_col': 2})
        membership.replace_column_cluster('C1', 1, 2)
        self.assertEquals([3, 2], membership.clusters_for_column('C1'))
        self.assertEquals(['C1', 'C2'], membership.columns_for_cluster(2))

    def test_replace_row_cluster(self):
        membership = memb.ClusterMembership(
            row_is_member_of={'R1': [1, 3], 'R2': [2, 3]},
            column_is_member_of={},
            config_params={'memb.clusters_per_row': 2})
        membership.replace_row_cluster('R1', 1, 2)
        self.assertEquals([3, 2], membership.clusters_for_row('R1'))
        self.assertEquals(['R1', 'R2'], membership.rows_for_cluster(2))

    def test_std_fuzzy_coefficient(self):
        result = [memb.std_fuzzy_coefficient(value, 10) for value in range(1, 11)]
        self.assertEquals([0.5685727544772026, 0.43416814526581843,
                           0.3345987618184194, 0.2608359483385415,
                           0.20619111210390084, 0.1657092217551106,
                           0.13571949977708733, 0.11350256730258876,
                           0.09704385891782485, 0.08485094785750477], result)

    def test_old_fuzzy_coefficient(self):
        result = [memb.old_fuzzy_coefficient(value, 2) for value in range(1, 3)]
        self.assertEquals([0.10150146242745953, 0.013736729166550634], result)


    def test_get_best_clusters(self):
        matrix = dm.DataMatrix(3, 4, values=[[1.0, 2.0, 3.0, 4.0],
                                          [4.0, 5.0, 6.0, 7.0],
                                          [7.0, 8.0, 9.0, 10.0]])
        result = memb.get_best_clusters(matrix, 2)
        self.assertTrue(len(result['Row 0']) == 2 and
                        3 in result['Row 0'] and 4 in result['Row 0'])
        self.assertTrue(len(result['Row 1']) == 2 and
                        3 in result['Row 1'] and 4 in result['Row 1'])
        self.assertTrue(len(result['Row 2']) == 2 and
                        3 in result['Row 2'] and 4 in result['Row 2'])

    def test_reseed_small_row_clusters(self):
        row_names = ['R1', 'R2', 'R3', 'R4', 'R5']
        membership = memb.ClusterMembership(
            row_is_member_of={'R1': [1, 3], 'R2': [2, 3]},
            column_is_member_of={},
            config_params={'memb.clusters_per_row': 2,
                           'memb.min_cluster_rows_allowed': 3,
                           'memb.max_cluster_rows_allowed': 70,
                           'num_clusters': 4})
        membership.reseed_small_row_clusters(row_names)
        for cluster in xrange(1, 5):
            self.assertTrue(membership.num_row_members(cluster) > 0)

    def test_reseed_small_column_clusters(self):
        col_names = ['C1', 'C2', 'C3']
        membership = memb.ClusterMembership(
            row_is_member_of={'R1': [1, 3], 'R2': [2, 3]},
            column_is_member_of={'C1': [1], 'C3': [2]},
            config_params={'memb.clusters_per_row': 2,
                           'memb.clusters_per_col': 2,
                           'num_clusters': 4})
        membership.reseed_small_column_clusters(col_names)
        for cluster in xrange(1, 5):
            self.assertTrue(membership.num_column_members(cluster) > 0)
        
