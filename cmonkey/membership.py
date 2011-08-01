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
            data_matrix, self.__row_membership, num_clusters,
            num_clusters_per_column)

        print "COLUMN MEMBERS: "
        for row in self.__column_membership:
            print row

    def cluster_for_row(self, row_index, column_index):
        """returns row membership value at the specified position"""
        return self.__row_membership[row_index][column_index]


def seed_column_members(data_matrix, row_membership, num_clusters,
                        num_clusters_per_column):
    """default column membership seeder"""
    num_rows = data_matrix.num_rows()
    num_cols = data_matrix.num_columns()
    print "# rows: %d, # clusters: %d" % (num_rows, num_clusters)
    # create a submatrix for each cluster
    column_scores = []
    for cluster_num in range(1, num_clusters + 1):
        current_cluster_rows = []
        for row_index in range(num_rows):
            if row_membership[row_index][0] == cluster_num:
                current_cluster_rows.append(row_index)
        submatrix = data_matrix.submatrix(current_cluster_rows,
                                          range(data_matrix.num_columns()))
        scores = -compute_column_scores(submatrix.values())
        print ("%d: " % cluster_num), scores
        column_scores.append(scores)

    column_members = []
    for column_index in range(num_cols):
        scores_to_order = []
        for row_index in range(num_clusters):
            scores_to_order.append(column_scores[row_index][column_index])
        column_members.append(order(scores_to_order)[:num_clusters_per_column])
    return column_members


def order(alist):
    """a weird R function that gives each item's position in the original list
    if you enumerate each item in a sorted list"""
    sorted_list = sorted(alist, reverse=True)
    return [(alist.index(item)) + 1 for item in sorted_list]


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
    return column_means(matrix_minus_colmeans_squared) / var_norm
