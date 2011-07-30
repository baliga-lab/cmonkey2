"""membership.py - cMonkey cluster membership management

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import numpy
from util import column_means


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

    def __init__(self, data_matrix,
                 num_clusters,
                 num_clusters_per_row,
                 num_clusters_per_column,
                 seed_row_memberships,
                 seed_column_memberships):
        """create instance of ClusterMembership"""
        num_rows = data_matrix.num_rows()
        self.__row_membership = [[0 for _ in range(num_clusters_per_row)]
                                 for _ in range(num_rows)]
        seed_row_memberships(self.__row_membership, data_matrix)
        self.__column_membership = seed_column_memberships(
            data_matrix, self.__row_membership, num_clusters)

    def cluster_for_row(self, row_index, column_index):
        """returns row membership value at the specified position"""
        return self.__row_membership[row_index][column_index]


def seed_column_members(data_matrix, row_membership, num_clusters):
    """default column membership seeder"""
    num_rows = data_matrix.num_rows()
    print "# rows: %d, # clusters: %d" % (num_rows, num_clusters)
    # create a submatrix for each cluster
    for cluster_num in range(1, num_clusters + 1):
        current_cluster_rows = []
        for row_index in range(num_rows):
            if row_membership[row_index][0] == cluster_num:
                current_cluster_rows.append(row_index)
        submatrix = data_matrix.submatrix(current_cluster_rows,
                                          range(data_matrix.num_columns()))
        column_scores = compute_column_scores(submatrix.values())
        print "COLUMN SCORES, CLUSTER %d" % cluster_num
        print column_scores


def compute_column_scores(matrix):
    """For a given matrix, compute the column scores.
    This is used to compute the column scores of the sub matrices that
    were determined by the pre-seeding, so typically, matrix is a
    submatrix of the input matrix that contains only the rows that
    belong to a certain cluster.

    This function normalizes diff^2 by the mean expression level, similar
    to "Index of Dispersion", see
    http://en.wikipedia.org/wiki/Index_of_dispersion
    for details
    """
    colmeans = column_means(matrix)
    matrix_minus_colmeans = []
    # subtract each row by the column means
    for row in matrix:
        new_row = []
        matrix_minus_colmeans.append(new_row)
        for col_index in range(len(row)):
            new_row.append(row[col_index] - colmeans[col_index])
    matrix_minus_colmeans_squared = numpy.square(matrix_minus_colmeans)
    var_norm = numpy.abs(colmeans) + 0.01
    return -(column_means(matrix_minus_colmeans_squared) / var_norm)
