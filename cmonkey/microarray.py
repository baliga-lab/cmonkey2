"""microarray.py - cMonkey microarray related processing

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import numpy
from util import column_means, row_means, quantile
from datamatrix import DataMatrix


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
    def __init__(self, row_is_member_of, column_is_member_of):
        """creates an instance of ClusterMembership"""

        def create_cluster_to_names_map(name_to_cluster_map):
            """from a name->cluster-list dictionary, create a cluster->name
            list dictionary"""
            result = {}
            for name, clusters in name_to_cluster_map.items():
                for cluster in clusters:
                    if cluster not in result:
                        result[cluster] = []
                    result[cluster].append(name)
            return result

        self.__row_is_member_of = row_is_member_of
        self.__column_is_member_of = column_is_member_of
        self.__cluster_row_members = create_cluster_to_names_map(
            self.__row_is_member_of)
        self.__cluster_column_members = create_cluster_to_names_map(
            self.__column_is_member_of)

    @classmethod
    def create(cls, data_matrix,
               num_clusters,
               num_clusters_per_row,
               num_clusters_per_column,
               seed_row_memberships,
               seed_column_memberships):
        """create instance of ClusterMembership using
        the provided seeding algorithms"""
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

        # using the seeding functions, build the initial membership
        # dictionaries
        num_rows = data_matrix.num_rows()
        row_membership = [[0 for _ in range(num_clusters_per_row)]
                          for _ in range(num_rows)]
        seed_row_memberships(row_membership, data_matrix)
        column_membership = seed_column_memberships(
            data_matrix, row_membership, num_clusters, num_clusters_per_column)
        row_is_member_of = make_member_map(row_membership,
                                           data_matrix.row_names())
        col_is_member_of = make_member_map(column_membership,
                                           data_matrix.column_names())
        return ClusterMembership(row_is_member_of, col_is_member_of)

    def clusters_for_row(self, row_name):
        """determine the clusters for the specified row"""
        return self.__row_is_member_of[row_name]

    def clusters_for_column(self, column_name):
        """determine the clusters for the specified column"""
        return self.__column_is_member_of[column_name]

    def rows_for_cluster(self, cluster):
        """determine the rows that belong to a cluster"""
        return self.__cluster_row_members[cluster]

    def columns_for_cluster(self, cluster):
        """determine the rows that belong to a cluster"""
        return self.__cluster_column_members[cluster]

    def is_row_member_of(self, row_name, cluster):
        """determines whether a certain row is member of a cluster"""
        return row_name in self.rows_for_cluster(cluster)


def seed_column_members(data_matrix, row_membership, num_clusters,
                        num_clusters_per_column):
    """Default column membership seeder
    In case of multiple input ratio matrices, we assume that these
    matrices have been combined into data_matrix"""
    num_rows = data_matrix.num_rows()
    num_cols = data_matrix.num_columns()
    # create a submatrix for each cluster
    column_scores = []
    for cluster_num in range(1, num_clusters + 1):
        current_cluster_rows = []
        for row_index in range(num_rows):
            if row_membership[row_index][0] == cluster_num:
                current_cluster_rows.append(data_matrix.row_name(row_index))
        submatrix = data_matrix.submatrix_by_name(
            row_names=current_cluster_rows)
        scores = (-compute_column_scores(submatrix))[0]
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
    The result is a DataMatrix with one row containing all the
    column scores

    This function normalizes diff^2 by the mean expression level, similar
    to "Index of Dispersion", see
    http://en.wikipedia.org/wiki/Index_of_dispersion
    for details
    """
    colmeans = matrix.column_means().values()[0]
    matrix_minus_colmeans_squared = __subtract_and_square(matrix, colmeans)
    var_norm = numpy.abs(colmeans) + 0.01
    result = column_means(matrix_minus_colmeans_squared) / var_norm
    return DataMatrix(1, matrix.num_columns(), ['Col. Scores'],
                      matrix.column_names(), [result])


def __subtract_and_square(matrix, vector):
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


def compute_row_scores(membership, matrix, num_clusters):
    """for each cluster 1, 2, .. num_clusters compute the row scores
    for the each row name in the input name matrix"""
    clusters = range(1, num_clusters + 1)
    cluster_row_scores = __compute_row_scores_for_clusters(
        membership, matrix, clusters)
    print "# rows: %d # cols: %d" % (cluster_row_scores[0].num_rows(),
                                     cluster_row_scores[0].num_columns())
    cluster_row_scores = __replace_non_numeric_values(cluster_row_scores,
                                                      membership,
                                                      matrix, clusters)

    # rearrange result into a DataMatrix, where rows are indexed by gene
    # and columns represent clusters
    result = DataMatrix(matrix.num_rows(), num_clusters,
                        row_names=matrix.row_names())
    for cluster in range(num_clusters):
        row_scores = cluster_row_scores[cluster]
        for row_index in range(matrix.num_rows()):
            result[row_index][cluster] = row_scores[0][row_index]
    print result.sorted_by_row_name()
    #for cluster in range(1, num_clusters + 1):
    #    print "CLUSTER %d" % cluster
    #    print cluster_row_scores


def __compute_row_scores_for_clusters(membership, matrix, clusters):
    """compute the pure row scores for the specified clusters
    without nowmalization"""
    result = []
    for cluster in clusters:
        sm1 = matrix.submatrix_by_name(
            row_names=membership.rows_for_cluster(cluster),
            column_names=membership.columns_for_cluster(cluster))
        if sm1.num_columns() > 1:
            matrix_filtered = matrix.submatrix_by_name(
                column_names=membership.columns_for_cluster(cluster))
            row_scores_for_cluster = __compute_row_scores_for_submatrix(
                matrix_filtered, sm1)
            result.append(row_scores_for_cluster)
        else:
            result.append(None)
    return result


def __compute_row_scores_for_submatrix(datamatrix, submatrix):
    """For a given matrix, compute the row scores. The second submatrix is
    used to calculate the column means on and should be derived from
    datamatrix filtered by the row names and column names of a specific
    cluster.
    datamatrix should be filtered by the columns of a specific cluster in
    order for the column means to be applied properly.
    The result is a DataMatrix with one row containing all the row scores"""
    colmeans = submatrix.column_means().values()[0]
    matrix_minus_colmeans_squared = __subtract_and_square(datamatrix, colmeans)
    scores = numpy.log(row_means(matrix_minus_colmeans_squared) + 1e-99)
    return DataMatrix(1, datamatrix.num_rows(),
                      row_names=['Row Scores'],
                      col_names=datamatrix.row_names(), values=[scores])


def __replace_non_numeric_values(cluster_row_scores, membership, matrix,
                                 clusters):
    """perform adjustments for NaN or inf values"""
    qvalue = __quantile_normalize_scores(cluster_row_scores, membership,
                                         clusters)
    result = []
    for row_scores in cluster_row_scores:
        if not row_scores:
            result.append(DataMatrix(1, matrix.num_rows(),
                                     col_names=matrix.row_names(),
                                     init_value=qvalue))
        else:
            for row_index in range(row_scores.num_rows()):
                for col_index in range(row_scores.num_columns()):
                    if not numpy.isfinite(row_scores[row_index][col_index]):
                        row_scores[row_index][col_index] = qvalue
            result.append(row_scores)
    return result


def __quantile_normalize_scores(cluster_row_scores, membership, clusters):
    """quantile normalize the row scores in cluster_row_scores
    that are not NaN or +/-Inf and are in a row cluster membership
    """
    values_for_quantile = []
    for cluster in clusters:
        row_scores_for_cluster = cluster_row_scores[cluster - 1]

        if row_scores_for_cluster != None:
            row_scores_names = row_scores_for_cluster.column_names()
            for col_index in range(row_scores_for_cluster.num_columns()):
                score = row_scores_for_cluster[0][col_index]
                gene_name = row_scores_names[col_index]
                if (numpy.isfinite(score)
                    and membership.is_row_member_of(gene_name, cluster)):
                    values_for_quantile.append(score)
    return quantile(values_for_quantile, 0.95)

__all__ = ['ClusterMembership', 'compute_row_scores', 'compute_column_scores',
           'seed_column_members']
