"""membership.py - cMonkey cluster membership management

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


class ClusterMembership:
    """A class to store row and column memberships of an input matrix
    1. Row memberships are stored as a matrix where each row represents a gene
    in the input and num_clusters_per_row columns.
    2. Column memberships are stored as a matrix of |input genes| rows and
       num_clusters_per_column columns columns
    When creating a ClusterMembership, a row membership seed function is
    called to generate the first column of cluster memberships
    This function is of the type
    DataMatrix => [int]
    A column seed function is called after this, which generates the
    entire column membership matrix.
    """

    def __init__(self, data_matrix, num_clusters_per_row,
                 num_clusters_per_column, seed_row_memberships,
                 seed_column_memberships):
        """create instance of ClusterMembership"""
        num_rows = data_matrix.num_rows()
        self.__row_membership = [[0.0 for _ in range(num_clusters_per_row)]
                                 for _ in range(num_rows)]
        seed_row_memberships(self.__row_membership, data_matrix)
        self.__column_membership = seed_column_memberships(
            self.__row_membership, data_matrix, num_clusters_per_column)

    def cluster_for_row(self, row_index, column_index):
        """returns row membership value at the specified position"""
        return self.__row_membership[row_index][column_index]
