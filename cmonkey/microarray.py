# vi: sw=4 ts=4 et:
"""microarray.py - cMonkey microarray related processing
This module captures the microarray-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import numpy as np
import logging
import cmonkey.datamatrix as dm
import cmonkey.util as util
import cmonkey.scoring as scoring

try:
    xrange
except NameError:
    xrange = range


def seed_column_members(data_matrix, row_membership, num_clusters,
                        num_clusters_per_column):
    """Default column membership seeder ('best')
    In case of multiple input ratio matrices, we assume that these
    matrices have been combined into data_matrix"""
    num_rows = data_matrix.num_rows
    num_cols = data_matrix.num_columns
    # create a submatrix for each cluster
    cscores = np.zeros([data_matrix.num_columns, num_clusters])
    for cluster_num in xrange(1, num_clusters + 1):
        current_cluster_rows = []
        for row_index in xrange(num_rows):
            if row_membership[row_index][0] == cluster_num:
                current_cluster_rows.append(data_matrix.row_names[row_index])
        submatrix = data_matrix.submatrix_by_name(
            row_names=current_cluster_rows)
        _, scores = scoring.compute_column_scores_submatrix(submatrix)
        cscores.T[cluster_num - 1] = -scores

    start_time = util.current_millis()
    column_members = [util.rorder(cscores[i], num_clusters_per_column)
                      for i in xrange(num_cols)]
    elapsed = util.current_millis() - start_time
    logging.debug("seed column members in %f s.", elapsed % 1000.0)
    return column_members


def compute_row_scores(membership, matrix, num_clusters, config_params):
    """for each cluster 1, 2, .. num_clusters compute the row scores
    for the each row name in the input name matrix"""
    start_time = util.current_millis()
    cluster_row_scores = __compute_row_scores_for_clusters(
        membership, matrix, num_clusters, config_params)
    # TODO: replace the nan/inf-Values with the quantile-thingy in the R-version

    logging.debug("__compute_row_scores_for_clusters() in %f s.",
                  (util.current_millis() - start_time) / 1000.0)

    # rearrange result into a DataMatrix, where rows are indexed by gene
    # and columns represent clusters
    start_time = util.current_millis()
    values = np.zeros((matrix.num_rows, num_clusters))

    # note that cluster is 0 based on a matrix
    for cluster in xrange(num_clusters):
        row_scores = cluster_row_scores[cluster]
        values[:, cluster] = row_scores
    result = dm.DataMatrix(matrix.num_rows, num_clusters,
                           row_names=matrix.row_names,
                           values=values)
    logging.debug("made result matrix in %f s.",
                  (util.current_millis() - start_time) / 1000.0)
    return result

ROW_SCORE_MATRIX = None
ROW_SCORE_MEMBERSHIP = None


def __compute_row_scores_for_clusters(membership, matrix, num_clusters,
                                      config_params):
    """compute the pure row scores for the specified clusters
    without nowmalization"""
    # note that we set the data into globals before we fork it off
    # to save memory and pickling time
    global ROW_SCORE_MATRIX, ROW_SCORE_MEMBERSHIP
    ROW_SCORE_MATRIX = matrix
    ROW_SCORE_MEMBERSHIP = membership

    if config_params['multiprocessing']:
        with util.get_mp_pool(config_params) as pool:
            result = pool.map(compute_row_scores_for_cluster, xrange(1, num_clusters + 1))
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

    if sm1.num_columns > 1:
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
    rm = util.row_means(np.square(matrix.values - util.column_means(submatrix.values)))
    # we clip the values to make sure the argument to log will be
    # sufficiently above 0 to avoid errors
    return np.log(np.clip(rm, 1e-20, 1000.0) + 1e-99)


class RowScoringFunction(scoring.ScoringFunctionBase):
    """Scoring algorithm for microarray data based on genes"""

    def __init__(self, function_id, cmrun):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, function_id, cmrun)
        self.run_log = scoring.RunLog(function_id, cmrun.dbsession(),
                                      cmrun.config_params)

    def do_compute(self, iteration_result, ref_matrix=None):
        """the row scoring function"""
        return compute_row_scores(self.membership,
                                  self.ratios,
                                  self.num_clusters(),
                                  self.config_params)

    def run_logs(self):
        """return the run logs"""
        return [self.run_log]

__all__ = ['compute_row_scores', 'seed_column_members']
