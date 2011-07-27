"""membership_test.py - unit test module for membership module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from membership import ClusterMembership


class MockDataMatrix:

    def __init__(self, num_rows):
        self.__num_rows = num_rows

    def num_rows(self):
        return self.__num_rows


class MockSeedRowMemberships:

    def __init__(self):
        self.was_called = False

    def __call__(self, row_membership, data_matrix):
        self.was_called = True


class MockSeedColumnMemberships:

    def __init__(self):
        self.was_called = False

    def __call__(self, row_membership, data_matrix, num_clusters_per_column):
        self.was_called = True
        return None


class ClusterMembershipTest(unittest.TestCase):
    """Test class for ClusterMembership"""

    def test_create_membership(self):
        """test creating a ClusterMembership object"""
        datamatrix = MockDataMatrix(3)
        seed_row_memberships = MockSeedRowMemberships()
        seed_col_memberships = MockSeedColumnMemberships()
        membership = ClusterMembership(datamatrix, 2, 2,
                                       seed_row_memberships,
                                       seed_col_memberships)
        self.assertTrue(seed_row_memberships.was_called)
        self.assertTrue(seed_col_memberships.was_called)
        self.assertEquals(0.0, membership.cluster_for_row(0, 1))
