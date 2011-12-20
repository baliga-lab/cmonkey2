"""microarray.py - cMonkey microarray related processing
This module captures the microarray-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import numpy as np
import logging
import datamatrix as dm
import util
import scoring
import multiprocessing as mp


def seed_column_members(data_matrix, row_membership, num_clusters,
                        num_clusters_per_column):
    """Default column membership seeder ('best')
    In case of multiple input ratio matrices, we assume that these
    matrices have been combined into data_matrix"""
    num_rows = data_matrix.num_rows()
    num_cols = data_matrix.num_columns()
    # create a submatrix for each cluster
    column_scores = []
    for cluster_num in xrange(1, num_clusters + 1):
        current_cluster_rows = []
        for row_index in xrange(num_rows):
            if row_membership[row_index][0] == cluster_num:
                current_cluster_rows.append(data_matrix.row_name(row_index))
        submatrix = data_matrix.submatrix_by_name(
            row_names=current_cluster_rows)
        scores = (-compute_column_scores_submatrix(submatrix))[0]
        column_scores.append(scores)

    column_members = []
    for column_index in xrange(num_cols):
        scores_to_order = []
        for row_index in xrange(num_clusters):
            scores_to_order.append(column_scores[row_index][column_index])
        column_members.append(order(scores_to_order)[:num_clusters_per_column])
    return column_members


def order(alist):
    """a weird R function that gives each item's position in the original list
    if you enumerate each item in a sorted list"""
    return [(alist.index(item)) + 1 for item in sorted(alist, reverse=True)]


def compute_column_scores(membership, matrix, num_clusters):
    """Computes the column scores for the specified number of clusters"""

    def compute_substitution(cluster_column_scores):
        """calculate substitution value for missing column scores"""
        membership_values = []
        for cluster in xrange(1, num_clusters + 1):
            columns = membership.columns_for_cluster(cluster)
            column_scores = cluster_column_scores[cluster - 1]
            if column_scores != None:
                for row in xrange(column_scores.num_rows()):
                    for col in xrange(column_scores.num_columns()):
                        if column_scores.column_name(col) in columns:
                            membership_values.append(column_scores[row][col])
        return util.quantile(membership_values, 0.95)

    cluster_column_scores = []
    null_scores_found = False
    for cluster in xrange(1, num_clusters + 1):
        submatrix = matrix.submatrix_by_name(
            row_names=membership.rows_for_cluster(cluster))
        if submatrix.num_rows() > 1:
            cluster_column_scores.append(compute_column_scores_submatrix(
                    submatrix))
        else:
            cluster_column_scores.append(None)
            null_scores_found = True

    if null_scores_found:
        substitution = compute_substitution(cluster_column_scores)

    # Convert scores into a matrix that have the clusters as columns
    # and conditions in the rows
    result = dm.DataMatrix(matrix.num_columns(), num_clusters,
                           row_names=matrix.column_names())
    for cluster in xrange(num_clusters):
        column_scores = cluster_column_scores[cluster]
        for row_index in xrange(matrix.num_columns()):
            if column_scores == None:
                result[row_index][cluster] = substitution
            else:
                result[row_index][cluster] = column_scores[0][row_index]
    return result


def compute_column_scores_submatrix(matrix):
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
    matrix_minus_colmeans_squared = subtract_and_square(matrix, colmeans)
    var_norm = np.abs(colmeans) + 0.01
    result = util.column_means(matrix_minus_colmeans_squared) / var_norm
    return dm.DataMatrix(1, matrix.num_columns(), ['Col. Scores'],
                         matrix.column_names(), [result])


def subtract_and_square(matrix, vector):
    """reusable function to subtract a vector from each row of
    the input matrix and square the values in the result matrix"""
    return np.square(matrix.values() - vector)


def compute_row_scores(membership, matrix, num_clusters,
                       use_multiprocessing):
    """for each cluster 1, 2, .. num_clusters compute the row scores
    for the each row name in the input name matrix"""
    #clusters = xrange(1, num_clusters + 1)
    start_time = util.current_millis()
    cluster_row_scores = __compute_row_scores_for_clusters(
        membership, matrix, num_clusters, use_multiprocessing)
    logging.info("__compute_row_scores_for_clusters() in %f s.",
                 (util.current_millis() - start_time) / 1000.0)
    start_time = util.current_millis()
    cluster_row_scores = __replace_non_numeric_values(cluster_row_scores,
                                                      membership,
                                                      matrix, num_clusters)
    logging.info("__replace_non_numeric_values() in %f s.",
                 (util.current_millis() - start_time) / 1000.0)

    # rearrange result into a DataMatrix, where rows are indexed by gene
    # and columns represent clusters
    start_time = util.current_millis()
    values = np.zeros((matrix.num_rows(), num_clusters))

    # note that cluster is 0 based on a matrix
    for cluster in xrange(num_clusters):
        row_scores = cluster_row_scores[cluster]
        values[:, cluster] = row_scores
    result = dm.DataMatrix(matrix.num_rows(), num_clusters,
                           row_names=matrix.row_names(),
                           values=values)
    logging.info("made result matrix in %f s.",
                 (util.current_millis() - start_time) / 1000.0)

    return result.sorted_by_row_name()

ROW_SCORE_MATRIX = None
ROW_SCORE_MEMBERSHIP = None


def __compute_row_scores_for_clusters(membership, matrix, num_clusters,
                                      use_multiprocessing):
    """compute the pure row scores for the specified clusters
    without nowmalization"""
    # note that we set the data into globals before we fork it off
    # to save memory and pickling time
    global ROW_SCORE_MATRIX, ROW_SCORE_MEMBERSHIP
    ROW_SCORE_MATRIX = matrix
    ROW_SCORE_MEMBERSHIP = membership

    if use_multiprocessing:
        pool = mp.Pool()
        result = pool.map(compute_row_scores_for_cluster,
                          [cluster for cluster in xrange(1, num_clusters + 1)])
        pool.close()
    else:
        result = []
        for cluster in range(1, num_clusters + 1):
            result.append(compute_row_scores_for_cluster(cluster))
    return result


def compute_row_scores_for_cluster(cluster):
    """This function computes the row score for a cluster"""
    global ROW_SCORE_MATRIX, ROW_SCORE_MEMBERSHIP
    membership = ROW_SCORE_MEMBERSHIP
    matrix = ROW_SCORE_MATRIX

    rnames = membership.rows_for_cluster(cluster)
    cnames = membership.columns_for_cluster(cluster)
    sm1 = matrix.submatrix_by_name(row_names=rnames, column_names=cnames)

    if sm1.num_columns() > 1:
        matrix_filtered = matrix.submatrix_by_name(column_names=cnames)
        row_scores_for_cluster = __compute_row_scores_for_submatrix(
            matrix_filtered, sm1)
        return row_scores_for_cluster
    else:
        return None


def __compute_row_scores_for_submatrix(matrix, submatrix):
    """For a given matrix, compute the row scores. The second submatrix is
    used to calculate the column means on and should be derived from
    datamatrix filtered by the row names and column names of a specific
    cluster.
    matrix should be filtered by the columns of a specific cluster in
    order for the column means to be applied properly.
    The result is a DataMatrix with one row containing all the row scores"""
    colmeans = submatrix.column_means().values()[0]
    return np.log(util.row_means(subtract_and_square(matrix, colmeans)) + 1e-99)


def __replace_non_numeric_values(cluster_row_scores, membership, matrix,
                                 num_clusters):
    """perform adjustments for NaN or inf values"""
    qvalue = __quantile_normalize_scores(cluster_row_scores,
                                         matrix.row_names(),
                                         membership,
                                         num_clusters)
    result = []
    for row_scores in cluster_row_scores:
        if row_scores == None:
            """no scores available, use the quantile normalized score"""
            row_scores = np.zeros(matrix.num_rows())
            row_scores.fill(qvalue)
        else:
            row_scores[np.isnan(row_scores)] = qvalue
        result.append(row_scores)

    return result


def __quantile_normalize_scores(cluster_row_scores,
                                row_names,
                                membership,
                                num_clusters):
    """quantile normalize the row scores in cluster_row_scores
    that are not NaN or +/-Inf and are in a row cluster membership
    """
    values_for_quantile = []
    for cluster in xrange(1, num_clusters + 1):
        row_scores_for_cluster = cluster_row_scores[cluster - 1]

        if row_scores_for_cluster != None:
            for row in xrange(len(row_scores_for_cluster)):
                score = row_scores_for_cluster[row]
                gene_name = row_names[row]
                if (np.isfinite(score)
                    and membership.is_row_member_of(gene_name, cluster)):
                    values_for_quantile.append(score)
    return util.quantile(values_for_quantile, 0.95)


class RowScoringFunction(scoring.ScoringFunctionBase):
    """Scoring algorithm for microarray data based on genes"""

    def __init__(self, membership, matrix, weight_func=None,
                 config_params=None):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, weight_func,
                                             config_params)

    def name(self):
        """returns the name of this scoring function"""
        return "Row"

    def compute(self, iteration, ref_matrix=None):
        """compute method, iteration is the 0-based iteration number"""
        start_time = util.current_millis()
        result = compute_row_scores(
            self.membership(),
            self.matrix(),
            self.num_clusters(),
            self.config_params[scoring.KEY_MULTIPROCESSING])

        elapsed = util.current_millis() - start_time
        logging.info("ROW SCORING TIME: %f s.", (elapsed / 1000.0))
        return result


class ColumnScoringFunction(scoring.ScoringFunctionBase):
    """Scoring algorithm for microarray data based on conditions.
    Note that the score does not correspond to the normal scoring
    function output format and can therefore not be combined in
    a generic way (the format is |condition x cluster|)"""

    def __init__(self, membership, matrix,
                 config_params=None):
        """create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, None, config_params)

    def name(self):
        """returns the name of this scoring function"""
        return "Column"

    def compute(self, iteration, ref_matrix=None):
        """compute method, iteration is the 0-based iteration number"""
        start_time = util.current_millis()
        result = compute_column_scores(self.membership(),
                                       self.matrix(),
                                       self.num_clusters())
        elapsed = util.current_millis() - start_time
        logging.info("COLUMN SCORING TIME: %f s.", (elapsed / 1000.0))
        return result

    def apply_weight(self, result, iteration):
        """applies the stored weight"""
        return result

__all__ = ['ClusterMembership', 'compute_row_scores', 'compute_column_scores',
           'seed_column_members']
