"""membership.py - cMonkey cluster membership management

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import numpy
from util import column_means, row_means


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
        # using the seeding functions, build the initial membership
        # dictionaries
        num_rows = data_matrix.num_rows()
        row_membership = [[0 for _ in range(num_clusters_per_row)]
                          for _ in range(num_rows)]
        seed_row_memberships(row_membership, data_matrix)
        column_membership = seed_column_memberships(data_matrix,
                                                    row_membership,
                                                    num_clusters,
                                                    num_clusters_per_column)
        self.__row_is_member_of = make_member_map(row_membership,
                                                  data_matrix.row_names())
        self.__col_is_member_of = make_member_map(column_membership,
                                                  data_matrix.column_names())


def make_member_map(membership, names):
    """using a membership array, build a dictionary representing
    the contained memberships for a name"""
    result = {}
    for row_index in range(len(names)):
        row = membership[row_index]
        result[names[row_index]] = sorted([row[col_index] for col_index
                                           in range(len(row))
                                           if row[col_index] > 0])
    return result


def seed_column_members(data_matrix, row_membership, num_clusters,
                        num_clusters_per_column):
    """default column membership seeder"""
    num_rows = data_matrix.num_rows()
    num_cols = data_matrix.num_columns()
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
    return [(alist.index(item)) + 1 for item in sorted(alist, reverse=True)]


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
    matrix_minus_colmeans_squared = compute_matrix_minus_vector_squared(
        matrix, colmeans)
    var_norm = numpy.abs(colmeans) + 0.01
    return column_means(matrix_minus_colmeans_squared) / var_norm


def compute_row_scores(matrix):
    """For a given matrix, compute the row scores"""
    colmeans = column_means(matrix)
    matrix_minus_colmeans_squared = compute_matrix_minus_vector_squared(
        matrix, colmeans)
    return numpy.log(row_means(matrix_minus_colmeans_squared) + 1e-99)


def compute_matrix_minus_vector_squared(matrix, vector):
    """reusable function to subtract a vector from each row of
    the input matrix and square the values in the result matrix"""
    result = []
    # subtract each row by the given vector
    for row in matrix:
        new_row = []
        result.append(new_row)
        for col_index in range(len(row)):
            new_row.append(row[col_index] - vector[col_index])
    return numpy.square(result)
