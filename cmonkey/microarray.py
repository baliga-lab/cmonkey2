# vi: sw=4 ts=4 et:
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
                current_cluster_rows.append(data_matrix.row_names[row_index])
        submatrix = data_matrix.submatrix_by_name(
            row_names=current_cluster_rows)
        scores = (-scoring.compute_column_scores_submatrix(submatrix)).values[0]
        column_scores.append(scores)

    column_members = []
    start_time = util.current_millis()
    for column_index in xrange(num_cols):
        scores_to_order = []
        for row_index in xrange(num_clusters):
            scores_to_order.append(column_scores[row_index][column_index])
        column_members.append(order(scores_to_order)[:num_clusters_per_column])
    elapsed = util.current_millis() - start_time
    logging.info("seed column members in %f s.", elapsed % 1000.0)
    return column_members

def order(alist):
    """a weird R function that gives each item's position in the original list
    if you enumerate each item in a sorted list"""
    return map(lambda x: alist.index(x) + 1, sorted(alist, reverse=True))


def compute_row_scores(membership, matrix, num_clusters,
                       use_multiprocessing):
    """for each cluster 1, 2, .. num_clusters compute the row scores
    for the each row name in the input name matrix"""
    start_time = util.current_millis()
    cluster_row_scores = __compute_row_scores_for_clusters(
        membership, matrix, num_clusters, use_multiprocessing)
    # TODO: replace the nan/inf-Values with the quantile-thingy in the R-version

    logging.info("__compute_row_scores_for_clusters() in %f s.",
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
                           row_names=matrix.row_names,
                           values=values)
    logging.info("made result matrix in %f s.",
                 (util.current_millis() - start_time) / 1000.0)

    result = result.sorted_by_row_name()
    result.fix_extreme_values()
    return result

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
        result = pool.map(compute_row_scores_for_cluster, xrange(1, num_clusters + 1))
        pool.close()
        pool.join()
    else:
        result = []
        for cluster in range(1, num_clusters + 1):
            result.append(compute_row_scores_for_cluster(cluster))
    # cleanup
    ROW_SCORE_MATRIX = None
    ROW_SCORE_MEMBERSHIP = None
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
    return np.log(
        util.row_means(np.square(matrix.values - submatrix.column_means())) + 1e-99)


"""
def __quantile_normalize_scores(cluster_row_scores,
                                row_names,
                                membership,
                                num_clusters):
    #quantile normalize the row scores in cluster_row_scores
    #that are not NaN or +/-Inf and are in a row cluster membership
    values_for_quantile = []
    for cluster in xrange(1, num_clusters + 1):
        row_scores_for_cluster = cluster_row_scores[cluster - 1]
        cluster_rows = membership.rows_for_cluster(cluster)
        if row_scores_for_cluster != None:
            for row in xrange(len(row_scores_for_cluster)):
                score = row_scores_for_cluster[row]
                gene_name = row_names[row]
                if np.isfinite(score) and (gene_name in cluster_rows):
                    values_for_quantile.append(score)
    return util.quantile(values_for_quantile, 0.95)
"""

class RowScoringFunction(scoring.ScoringFunctionBase):
    """Scoring algorithm for microarray data based on genes"""

    def __init__(self, membership, matrix, scaling_func=None,
                 run_in_iteration=scoring.schedule(1, 2),
                 config_params=None):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, scaling_func,
                                             run_in_iteration,
                                             config_params)
        self.cache_result = True
        self.run_log = scoring.RunLog("row_scoring", config_params)

    def name(self):
        """returns the name of this scoring function"""
        return "Row"

    def do_compute(self, iteration_result, ref_matrix=None):
        """the row scoring function"""
        return compute_row_scores(self.membership(),
                                  self.matrix(),
                                  self.num_clusters(),
                                  self.config_params[scoring.KEY_MULTIPROCESSING])

    def run_logs(self):
        """return the run logs"""
        return [self.run_log]

__all__ = ['compute_row_scores', 'seed_column_members']
