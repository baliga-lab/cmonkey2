"""microarray_test.py - unit test module for microarray module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from microarray import ClusterMembership, compute_column_scores
from microarray import seed_column_members
from datamatrix import DataMatrix
import numpy


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
        self.__num_rows = num_rows

    def num_rows(self):
        return self.__num_rows

    def num_columns(self):
        return 5

    def row_names(self):
        return [('GENE%d' % index) for index in range(self.__num_rows)]

    def column_names(self):
        return [('COND%d' % index) for index in range(self.num_columns())]

    def values(self):
        return MATRIX1

    def submatrix(self, rows, cols):
        return self


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


class ClusterMembershipTest(unittest.TestCase):
    """Test class for ClusterMembership"""

    def test_constructor(self):
        """tests the constructor"""
        membership = ClusterMembership({'R1': [1, 3], 'R2': [2, 3]},
                                       {'C1': [1, 2], 'C2': [2]})
        self.assertEquals([1, 3], membership.clusters_for_row('R1'))
        self.assertEquals([2, 3], membership.clusters_for_row('R2'))
        self.assertEquals([1, 2], membership.clusters_for_column('C1'))
        self.assertEquals([2], membership.clusters_for_column('C2'))
        self.assertEquals(['R1'], membership.rows_for_cluster(1))
        self.assertEquals(['R2'], membership.rows_for_cluster(2))
        self.assertEquals(['R1', 'R2'], membership.rows_for_cluster(3))
        self.assertEquals(['C1'], membership.columns_for_cluster(1))
        self.assertTrue(membership.is_row_member_of('R1', 1))
        self.assertTrue(membership.is_row_member_of('R2', 2))
        self.assertFalse(membership.is_row_member_of('R1', 2))

    def test_create_membership(self):
        """test creating a ClusterMembership object"""
        datamatrix = MockDataMatrix(3)
        seed_row_memberships = MockSeedRowMemberships()
        seed_col_memberships = MockSeedColumnMemberships()
        membership = ClusterMembership.create(datamatrix, 1, 2, 2,
                                              seed_row_memberships,
                                              seed_col_memberships)
        self.assertTrue(seed_row_memberships.was_called)
        self.assertTrue(seed_col_memberships.was_called)

    def test_compute_column_scores(self):
        """tests compute_column_scores"""
        matrix = DataMatrix(10, 5, ['R1', 'R2', 'R3', 'R4', 'R5', 'R6',
                                    'R7', 'R8', 'R9', 'R10'],
                            ['C1', 'C2', 'C3', 'C4', 'C5'],
                            MATRIX1)
        result = compute_column_scores(matrix)
        scores = result[0]
        self.assertEqual(5, len(scores))
        self.assertAlmostEqual(0.03085775, scores[0])
        self.assertAlmostEqual(0.05290099, scores[1])
        self.assertAlmostEqual(0.05277032, scores[2])
        self.assertAlmostEqual(0.00358045, scores[3])
        self.assertAlmostEqual(0.03948821, scores[4])

    def test_seed_column_members(self):
        """tests seed_column_members"""
        data_matrix = DataMatrix(3, 2, ["GENE1", "GENE2", "GENE3"],
                                 ['COL1', 'COL2'], [[1.0, 2.0], [2.0, 1.0],
                                                    [2.0, 1.0]])
        row_membership = [[1, 0], [2, 0], [1, 0]]
        column_members = seed_column_members(data_matrix, row_membership,
                                             2, 2)
        self.assertEquals([2, 1], column_members[0])
        self.assertEquals([2, 1], column_members[1])
