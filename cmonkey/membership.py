# vi: sw=4 ts=4 et:
"""membership.py - cMonkey cluster membership functionality
This module captures the microarray-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import math
import random
import logging
import sys
import numpy as np
import rpy2.robjects as robjects
from sqlalchemy import func

import cmonkey.datamatrix as dm
import cmonkey.util as util

# Python2/Python3 compatibility
try:
    import cPickle as pickle
except ImportError:
    import pickle

try:
    xrange
except NameError:
    xrange = range

import array
from collections import defaultdict


# Default values for membership creation
MAX_ADJUST_TRIES = 50

KEY_NUM_CLUSTERS = 'num_clusters'
KEY_CLUSTERS_PER_ROW = 'memb.clusters_per_row'
KEY_CLUSTERS_PER_COL = 'memb.clusters_per_col'
KEY_PROB_ROW_CHANGE = 'memb.prob_row_change'
KEY_PROB_COL_CHANGE = 'memb.prob_col_change'
KEY_MAX_CHANGES_PER_ROW = 'memb.max_changes_per_row'
KEY_MAX_CHANGES_PER_COL = 'memb.max_changes_per_col'
KEY_MIN_CLUSTER_ROWS_ALLOWED = 'memb.min_cluster_rows_allowed'
KEY_MAX_CLUSTER_ROWS_ALLOWED = 'memb.max_cluster_rows_allowed'

# These keys are for save points
KEY_ROW_IS_MEMBER_OF = 'memb.row_is_member_of'
KEY_COL_IS_MEMBER_OF = 'memb.col_is_member_of'


class OrigMembership:
    """This is an implementation of a membership data structure that more
    closely resembles the R original. It is much simpler than
    ClusterMembership, with a smaller memory footprint"""
    def __init__(self, row_names, col_names,
                 row_is_member_of, col_is_member_of,
                 config_params, row_indexes=None, col_indexes=None):
        """identical constructor to ClusterMembership"""
        self.__config_params = config_params

        # table with |genes| rows and the configured number of columns
        num_per_row = config_params['memb.clusters_per_row']
        num_per_col = config_params['memb.clusters_per_col']
        self.row_names = row_names
        self.col_names = col_names

        # reuse column indexes if they were provided
        if row_indexes:
            self.rowidx = row_indexes
        else:
            self.rowidx = {name: i for i, name in enumerate(self.row_names)}
        if col_indexes:
            self.colidx = col_indexes
        else:
            self.colidx = {name: i for i, name in enumerate(self.col_names)}

        """
        main representation:
        A matrix that has genes/conditions as their row, in the order of row_names/col_names.
        Each row in the matrix has num_per_X slots, each cell contains a cluster number.
        """
        self.row_membs = np.zeros((len(row_names), num_per_row), dtype='int32')
        self.col_membs = np.zeros((len(col_names), num_per_col), dtype='int32')

        for row, clusters in row_is_member_of.items():
            tmp = row_is_member_of[row][:num_per_row]
            for i in range(len(tmp)):
                self.row_membs[self.rowidx[row]][i] = tmp[i]

        for col, clusters in col_is_member_of.items():
            tmp = col_is_member_of[col][:num_per_col]
            for i in range(len(tmp)):
                self.col_membs[self.colidx[col]][i] = tmp[i]

    def write_column_members(self, filename):
        """Mostly for debugging, write out the current column membership state into a TSV file"""
        with open(filename, 'w') as outfile:
            colnums = range(1, len(self.col_membs[0]) + 1)
            outfile.write('\t'.join(map(lambda i: 'V%d' % i, colnums)))
            outfile.write('\n')
            for colname in sorted(self.col_names):
                idx = self.colidx[colname]
                row = self.col_membs[idx]
                outfile.write('%s\t' % colname)
                outfile.write('\t'.join(map(str, row)))
                outfile.write('\n')

    def write_row_members(self, filename):
        """Mostly for debugging, write out the current row membership state into a TSV file"""
        with open(filename, 'w') as outfile:
            colnums = range(1, len(self.row_membs[0]) + 1)
            outfile.write('\t'.join(map(lambda i: 'V%d' % i, colnums)))
            outfile.write('\n')
            for rowname in sorted(self.row_names):
                idx = self.rowidx[rowname]
                row = self.row_membs[idx]
                outfile.write('%s\t' % rowname)
                outfile.write('\t'.join(map(str, row)))
                outfile.write('\n')

    def num_clusters(self):
        """returns the number of clusters"""
        return self.__config_params[KEY_NUM_CLUSTERS]

    def num_clusters_per_row(self):
        """returns the number of clusters per row"""
        return self.__config_params[KEY_CLUSTERS_PER_ROW]

    def num_clusters_per_column(self):
        """returns the number of clusters per row"""
        return self.__config_params[KEY_CLUSTERS_PER_COL]

    def probability_seeing_row_change(self):
        """returns the probability for seeing a row change"""
        return self.__config_params[KEY_PROB_ROW_CHANGE]

    def probability_seeing_col_change(self):
        """returns the probability for seeing a row change"""
        return self.__config_params[KEY_PROB_COL_CHANGE]

    def max_changes_per_row(self):
        """returns the maximum number of changes per row"""
        return self.__config_params[KEY_MAX_CHANGES_PER_ROW]

    def max_changes_per_col(self):
        """returns the maximum number of changes per column"""
        return self.__config_params[KEY_MAX_CHANGES_PER_COL]

    def min_cluster_rows_allowed(self):
        """returns the minimum number of rows that should be in a cluster"""
        return self.__config_params[KEY_MIN_CLUSTER_ROWS_ALLOWED]

    def max_cluster_rows_allowed(self):
        """returns the maximum number of rows that should be in a cluster"""
        return self.__config_params[KEY_MAX_CLUSTER_ROWS_ALLOWED]

    def min_cluster_columns_allowed(self):
        """returns the minimum number of columns that should be in a cluster"""
        return 0

    def clusters_for_row(self, row):
        """determine the clusters for the specified row"""
        c = self.row_membs[self.rowidx[row]]
        return set(c[c > 0])

    def num_clusters_for_row(self, row):
        """returns the number of clusters for the row"""
        return len(self.clusters_for_row(row))

    def clusters_for_column(self, column):
        """determine the clusters for the specified column"""
        c = self.col_membs[self.colidx[column]]
        return set(c[c > 0])

    def num_clusters_for_column(self, column):
        """returns the number of clusters for the column"""
        return len(self.clusters_for_column(column))

    def rows_for_cluster(self, cluster):
        idx = set(np.where(self.row_membs == cluster)[0])
        return {self.row_names[i] for i in idx}

    def columns_for_cluster(self, cluster):
        idx = set(np.where(self.col_membs == cluster)[0])
        return {self.col_names[i] for i in idx}

    def num_row_members(self, cluster):
        return len(self.rows_for_cluster(cluster))

    def num_column_members(self, cluster):
        return len(self.columns_for_cluster(cluster))

    def clusters_not_in_row(self, row, clusters):
        return [cluster for cluster in clusters
                if cluster not in self.clusters_for_row(row)]

    def clusters_not_in_column(self, col, clusters):
        return [cluster for cluster in clusters
                if cluster not in self.clusters_for_column(col)]

    def is_row_in_cluster(self, row, cluster):
        return cluster in self.clusters_for_row(row)

    def is_column_in_cluster(self, col, cluster):
        return cluster in self.clusters_for_column(col)

    def free_slots_for_row(self, row):
        return np.where(self.row_membs[self.rowidx[row]] == 0)[0]

    def free_slots_for_column(self, col):
        return np.where(self.col_membs[self.colidx[col]] == 0)[0]

    def add_cluster_to_row(self, row, cluster, force=False):
        rowidx = self.rowidx[row]
        free_slots = np.where(self.row_membs[rowidx] == 0)[0]
        if len(free_slots > 0):
            index = free_slots[0]
            self.row_membs[rowidx, index] = cluster
        elif not force:
            raise Exception(("add_cluster_to_row() - exceeded clusters/row " +
                             "limit for row: '%s'" % str(row)))
        else:
            tmp = np.zeros((self.row_membs.shape[0], self.row_membs.shape[1] + 1), dtype='int32')
            tmp[:, :-1] = self.row_membs
            self.row_membs = tmp
            self.row_membs[rowidx][-1] = cluster

    def add_cluster_to_column(self, col, cluster, force=False):
        colidx = self.colidx[col]
        free_slots = np.where(self.col_membs[colidx] == 0)[0]
        if len(free_slots) > 0:
            index = free_slots[0]
            self.col_membs[colidx, index] = cluster
        elif not force:
            raise Exception(("add_cluster_to_column() - exceeded clusters/col " +
                             "limit for column: '%s'" % str(col)))
        else:
            tmp = np.zeros((self.col_membs.shape[0], self.col_membs.shape[1] + 1), dtype='int32')
            tmp[:, :-1] = self.col_membs
            self.col_membs = tmp
            self.col_membs[colidx][-1] = cluster

    def replace_row_cluster(self, row, index, new):
        self.row_membs[self.rowidx[row], index] = new

    def replace_column_cluster(self, col, index, new):
        self.col_membs[self.colidx[col], index] = new

    def pickle_path(self):
        """returns the function-specific pickle-path"""
        return '%s/last_row_scores.pkl' % (self.__config_params['output_dir'])

    def update(self, matrix, row_scores, column_scores,
               num_iterations, iteration_result):
        """top-level update method"""
        start = util.current_millis()
        row_scores, column_scores = fuzzify(self, row_scores, column_scores,
                                            num_iterations, iteration_result,
                                            self.__config_params['add_fuzz'])
        elapsed = util.current_millis() - start
        logging.debug("fuzzify took %f s.", elapsed / 1000.0)

        # pickle the (potentially fuzzed) row scores to use them
        # in the post adjustment step. We only need to do that in the last
        # iteration
        iteration = iteration_result['iteration']
        if iteration == num_iterations:
            with open(self.pickle_path(), 'wb') as outfile:
                pickle.dump(row_scores, outfile)

        start = util.current_millis()
        rd_scores, cd_scores = get_density_scores(self, row_scores,
                                                  column_scores)
        elapsed = util.current_millis() - start
        logging.debug("GET_DENSITY_SCORES() took %f s.", elapsed / 1000.0)

        start = util.current_millis()
        compensate_size(self, matrix, rd_scores, cd_scores)
        elapsed = util.current_millis() - start
        logging.debug("COMPENSATE_SIZE() took %f s.", elapsed / 1000.0)

        start_time = util.current_millis()
        update_for_rows(self, rd_scores, self.__config_params['multiprocessing'])
        elapsed = util.current_millis() - start_time
        logging.debug("update_for rdscores finished in %f s.", elapsed / 1000.0)

        start_time = util.current_millis()
        update_for_cols(self, cd_scores, self.__config_params['multiprocessing'])
        elapsed = util.current_millis() - start_time
        logging.debug("update_for cdscores finished in %f s.", elapsed / 1000.0)


def create_membership(matrix, seed_row_memberships, seed_column_memberships,
                      config_params):
    """create instance of ClusterMembership using
    the provided seeding algorithms"""
    def make_member_map(membs, names):
        """build a map row->[clusters]"""
        result = {}
        for i in xrange(len(names)):
            result[names[i]] = [c for c in membs[i] if c > 0]
        return result

    # using the seeding functions, build the initial membership
    # dictionaries
    num_clusters_per_row = config_params[KEY_CLUSTERS_PER_ROW]
    num_clusters = config_params[KEY_NUM_CLUSTERS]
    num_clusters_per_col = config_params[KEY_CLUSTERS_PER_COL]

    num_rows = matrix.num_rows
    row_membership = [[0 for _ in xrange(num_clusters_per_row)]
                      for _ in xrange(num_rows)]
    seed_row_memberships(row_membership, matrix)
    column_membership = seed_column_memberships(matrix, row_membership,
                                                num_clusters, num_clusters_per_col)
    row_is_member_of = make_member_map(row_membership, matrix.row_names)
    col_is_member_of = make_member_map(column_membership, matrix.column_names)
    return OrigMembership(matrix.row_names, matrix.column_names,
                          row_is_member_of, col_is_member_of,
                          config_params, matrix.row_indexes, matrix.column_indexes)


def update_for_rows(membership, rd_scores, multiprocessing):
    """generically updating row memberships according to  rd_scores"""
    rownames = rd_scores.row_names
    # note: for rows, the original version sorts the best clusters by cluster number !!!
    best_clusters = get_best_clusters(rd_scores, membership.num_clusters_per_row(), True)
    max_changes = membership.max_changes_per_row()
    change_prob = membership.probability_seeing_row_change()

    for index in xrange(rd_scores.num_rows):
        row = rownames[index]
        clusters = best_clusters[row]

        if seeing_change(change_prob):
            for _ in range(max_changes):
                if len(clusters) > 0:
                    free_slots = membership.free_slots_for_row(row)
                    if len(free_slots) > 0:
                        take_cluster = clusters[free_slots[0]]
                        if take_cluster not in membership.clusters_for_row(row):
                            membership.add_cluster_to_row(row, take_cluster)
                    else:
                        replace_delta_row_member(membership, row, clusters, rd_scores)


def replace_delta_row_member(membership, row, rm, rd_scores):
    index = rd_scores.row_indexes_for([row])[0]
    rds_values = rd_scores.values
    # Since Python is 0-based, we adjust the clusters by -1 to access the
    # arrays. This a little confusing, so we need to pay attention to this
    # function
    curr_clusters = membership.row_membs[membership.rowidx[row]] - 1
    rm_clusters = np.array(rm)
    rm_clusters -= 1
    deltas = rds_values[index, rm_clusters] - rds_values[index, curr_clusters]

    # ignore the positions in curr_cluster that are also in rm_clusters
    # delta 0 is a non-replacement
    erase = [i for i, cluster in enumerate(curr_clusters) if cluster in rm_clusters]
    for i in erase:
        deltas[i] = 0
    if len(deltas[deltas != 0.0]) > 0:
        maxidx = deltas.argmax(axis=0)
        # Note: this means clusters can only be assigned to rows once, we don't
        # really have to check whether we have more than 1 of the same cluster
        if rm[maxidx] not in membership.clusters_for_row(row):
            membership.replace_row_cluster(row, maxidx, rm[maxidx])


def update_for_cols(membership, cd_scores, multiprocessing):
    """updating column memberships according to cd_scores"""
    global UPDATE_MEMBERSHIP

    colnames = cd_scores.row_names
    best_clusters = get_best_clusters(cd_scores, membership.num_clusters_per_column())
    max_changes = membership.max_changes_per_col()
    change_prob = membership.probability_seeing_col_change()

    for index in xrange(cd_scores.num_rows):
        col = colnames[index]
        clusters = best_clusters[col]
        if seeing_change(change_prob):
            for c in range(max_changes):
                if len(clusters) > 0:
                    free_slots = membership.free_slots_for_column(col)
                    if len(free_slots) > 0:
                        slot = free_slots[0]
                        """
                        the slot can be out of bounds for the clusters array when
                        the setting for clusters_per_row/clusters_per_col is too
                        large, in this case pick a spot inside the array to avoid
                        the exception"""
                        if slot > len(clusters) - 1:
                            slot = len(clusters) - 1
                        take_cluster = clusters[slot]
                        #print "ii = ", c, ", add cluster ", take_cluster, " at ", free_slots[0]
                        membership.add_cluster_to_column(col, take_cluster)
                    else:
                        col_clusters = membership.col_membs[membership.colidx[col]]
                        multi = util.which_multiple(col_clusters)
                        if len(multi) > 0:
                            # indexes of col_clusters that are in multiple
                            for i, cluster in enumerate(col_clusters):
                                if cluster in multi:
                                    #print "multiple in row: ", index, " ii: ", c, " col.change ", i, " -> ", clusters[i]
                                    membership.replace_column_cluster(col, i, clusters[i])
                                    break
                        else:
                            replace_delta_column_member(membership, col, clusters, cd_scores)


def replace_delta_column_member(membership, col, cm, cd_scores):
    index = cd_scores.row_indexes_for([col])[0]
    cds_values = cd_scores.values
    curr_clusters = membership.col_membs[membership.colidx[col]] - 1
    cm_clusters = np.array(cm)
    cm_clusters -= 1
    deltas = cds_values[index, cm_clusters] - cds_values[index, curr_clusters]

    if len(deltas[deltas != 0.0]) > 0:
        maxidx = deltas.argmax(axis=0)
        #print "replace_delta col.change ", maxidx, " -> ", cm[maxidx]
        # Note: columns allow multiple cluster assignment !!!
        membership.replace_column_cluster(col, maxidx, cm[maxidx])


def postadjust(membership, rowscores, cutoff=0.33, limit=100):
    """adjusting the cluster memberships after the main iterations have been done
    Returns true if the function changed the membership, false if not"""
    assign_list = []
    for cluster in range(1, membership.num_clusters() + 1):
        assign = adjust_cluster(membership, cluster, rowscores, cutoff, limit)
        assign_list.append(assign)

    for assign in assign_list:
        for row, cluster in assign.items():
            membership.add_cluster_to_row(row, cluster, force=True)


def adjust_cluster(membership, cluster, rowscores, cutoff, limit):
    """adjust a single cluster"""
    def max_row_in_column(matrix, column):
        """returns a pair of the maximum row index and score in the given matrix and column"""
        sm = matrix.submatrix_by_name(wh, [matrix.column_names[column]])
        sm_values = sm.values
        max_row = 0
        max_score = -sys.float_info.max
        for row in range(sm.num_rows):
            if sm_values[row][0] > max_score:
                max_score = sm_values[row, 0]
                max_row = row
        return sm.row_names[max_row]

    old_rows = membership.rows_for_cluster(cluster)
    not_in = [(i, row) for i, row in enumerate(rowscores.row_names)
              if row not in old_rows]
    threshold = rowscores.submatrix_by_name(old_rows,
                                            [rowscores.column_names[cluster - 1]]).quantile(cutoff)
    wh = []
    rs_values = rowscores.values
    for row, row_name in not_in:
        if rs_values[row, cluster - 1] < threshold:
            wh.append(row_name)
    if len(wh) == 0 or len(wh) > limit:
        return {}  # return unmodified row membership

    tries = 0
    result = {}
    while len(wh) > 0 and tries < MAX_ADJUST_TRIES:
        wh2 = max_row_in_column(rowscores, cluster - 1)
        result[wh2] = cluster
        wh.remove(wh2)
        tries += 1
    old_num = len(membership.rows_for_cluster(cluster))
    logging.debug("CLUSTER %d, # ROWS BEFORE: %d, AFTER: %d",
                  cluster, old_num, old_num + len(result))
    return result

######################################################################
### Helpers
######################################################################

# Parallelized updating of row and column membership changes
UPDATE_MEMBERSHIP = None


def seeing_change(prob):
    """returns true if the update is seeing the change"""
    return prob >= 1.0 or random.uniform(0.0, 1.0) <= prob


def get_best_clusters(scores, n, sort=False):
    """retrieve the n best scored clusters for the given row/column score matrix"""
    if sort:
        return {scores.row_names[row]: sorted(util.rorder(scores.row_values(row), n))
                for row in xrange(scores.num_rows)}
    else:
        return {scores.row_names[row]: util.rorder(scores.row_values(row), n)
                for row in xrange(scores.num_rows)}


def get_row_density_scores(membership, row_scores):
    """getting density scores improves small clusters"""
    num_clusters = membership.num_clusters()
    rscore_range = abs(row_scores.max() - row_scores.min())
    rowscore_bandwidth = max(rscore_range / 100.0, 0.001)
    rd_scores = dm.DataMatrix(row_scores.num_rows,
                              row_scores.num_columns,
                              row_scores.row_names,
                              row_scores.column_names)
    rds_values = rd_scores.values

    start_time = util.current_millis()
    for cluster in xrange(1, num_clusters + 1):
        # instead of assigning the rr_scores values per row, we can assign to the
        # transpose and let numpy do the assignment
        rds_values.T[cluster - 1] = get_rr_scores(membership, row_scores,
                                                  rowscore_bandwidth,
                                                  cluster)

    elapsed = util.current_millis() - start_time
    logging.debug("RR_SCORES IN %f s.", elapsed / 1000.0)
    return rd_scores


def get_col_density_scores(membership, col_scores):
    num_clusters = membership.num_clusters()
    cscore_range = abs(col_scores.max() - col_scores.min())
    colscore_bandwidth = max(cscore_range / 100.0, 0.001)
    cd_scores = dm.DataMatrix(col_scores.num_rows,
                              col_scores.num_columns,
                              col_scores.row_names,
                              col_scores.column_names)
    cds_values = cd_scores.values

    start_time = util.current_millis()
    for cluster in xrange(1, num_clusters + 1):
        # instead of assigning the cc_scores values per row, we can assign to the
        # transpose and let numpy do the assignment
        cds_values.T[cluster - 1] = get_cc_scores(membership, col_scores,
                                                  colscore_bandwidth,
                                                  cluster)

    elapsed = util.current_millis() - start_time
    logging.debug("CC_SCORES IN %f s.", elapsed / 1000.0)
    return cd_scores


def get_density_scores(membership, row_scores, col_scores):
    return (get_row_density_scores(membership, row_scores),
            get_col_density_scores(membership, col_scores))


def get_rr_scores(membership, rowscores, bandwidth, cluster):
    """calculate the density scores for the given row score values in the
    specified cluster"""
    def bwscale(value):
        """standard bandwidth scaling function for row scores"""
        return math.exp(-value / 10.0) * 10.0

    cluster_rows = membership.rows_for_cluster(cluster)
    cluster_columns = membership.columns_for_cluster(cluster)
    kscores = rowscores.column_values(cluster - 1)
    kscores_finite = kscores[np.isfinite(kscores)]

    if len(cluster_rows) == 0 or len(kscores_finite) == 0 or len(cluster_columns) == 0:
        num_rows = rowscores.num_rows
        return [(1.0 / num_rows) for _ in xrange(num_rows)]
    else:
        score_indexes = rowscores.row_indexes_for(cluster_rows)
        cluster_scores = [kscores[index] for index in score_indexes]
        cluster_bandwidth = bandwidth * bwscale(len(cluster_rows))
        return util.density(kscores, cluster_scores, cluster_bandwidth,
                            np.amin(kscores_finite) - 1,
                            np.amax(kscores_finite) + 1)


def get_cc_scores(membership, scores, bandwidth, cluster):
    """calculate the density scores for the given column score values in the
    specified cluster"""
    cluster_rows = membership.rows_for_cluster(cluster)
    cluster_columns = membership.columns_for_cluster(cluster)
    kscores = scores.column_values(cluster - 1)
    kscores_finite = kscores[np.isfinite(kscores)]

    if len(cluster_rows) == 0 or len(kscores_finite) == 0 or len(cluster_columns) <= 1:
        # This is a little weird, but is here to at least attempt to simulate
        # what the original cMonkey is doing
        num_rows = scores.num_rows
        return [(1.0 / num_rows) for _ in xrange(num_rows)]
    else:
        score_indexes = scores.row_indexes_for(cluster_columns)
        cluster_scores = [kscores[index] for index in score_indexes]
        return util.density(kscores, cluster_scores, bandwidth,
                            np.amin(kscores_finite) - 1,
                            np.amax(kscores_finite) + 1)


def compensate_size(membership, matrix, rd_scores, cd_scores):
    """size compensation function"""
    def compensate_dim_size(size, dimsize, clusters_per_dim, num_clusters):
        """compensate size for a dimension"""
        return math.exp(-float(size) / (float(dimsize) *
                                        float(clusters_per_dim) /
                                        float(num_clusters)))

    def compensate_row_size(size):
        """compensation function for row dimension"""
        return compensate_dim_size(size,
                                   matrix.num_rows,
                                   membership.num_clusters_per_row(),
                                   membership.num_clusters())

    def compensate_column_size(size):
        """compensation function for column dimension"""
        return compensate_dim_size(size,
                                   matrix.num_columns,
                                   membership.num_clusters_per_column(),
                                   membership.num_clusters())

    def compensate_rows(cluster):
        """compensate density scores for row dimension"""
        num_rowmembers = membership.num_row_members(cluster)
        rd_scores.multiply_column_by(
            cluster - 1,
            compensate_row_size(max(num_rowmembers,
                                    membership.min_cluster_rows_allowed())))

    def compensate_columns(cluster):
        """compensate density scores for column dimension"""
        num_colmembers = membership.num_column_members(cluster)
        cd_scores.multiply_column_by(
            cluster - 1,
            compensate_column_size(max(num_colmembers,
                                       matrix.num_columns / 10.0)))

    num_clusters = membership.num_clusters()
    for cluster in xrange(1, num_clusters + 1):
        compensate_rows(cluster)
        compensate_columns(cluster)


def std_fuzzy_coefficient(iteration, num_iterations):
    """standard fuzzy coefficient as defined in cMonkey"""
    return 0.7 * math.exp(-(float(iteration) /
                            (float(num_iterations) / 3.0))) + 0.05


def old_fuzzy_coefficient(iteration, num_iterations):
    """standard fuzzy coefficient as defined in cMonkey"""
    return 0.75 * math.exp(-iteration/(num_iterations/4.0))


def make_kmeans_row_seeder(num_clusters):
    """creates a row seeding function based on k-means"""

    def seed(row_membership, matrix):
        """uses k-means seeding to seed row membership"""
        flat_values = matrix.values.flatten()
        flat_values[np.isnan(flat_values)] = 0.0
        flat_values[np.isinf(flat_values)] = 0.0
        flat_values[np.isneginf(flat_values)] = 0.0
        matrix_values = robjects.r.matrix(
            robjects.FloatVector(flat_values), nrow=matrix.num_rows, byrow=True)
        kmeans = robjects.r['kmeans']
        kwargs = {'centers': num_clusters, 'iter.max': 20, 'nstart': 2}
        seeding = kmeans(matrix_values, **kwargs)[0]
        for row in xrange(len(seeding)):
            row_membership[row][0] = seeding[row]

    return seed


def make_file_seeder(filename, sep=' '):
    """uses a TSV file to seed row membership"""

    def seed(row_membership, matrix):
        """ignore matrix parameter"""
        row_map = {name: idx
                   for idx, name in enumerate(matrix.row_names)}

        with open(filename) as infile:
            header = infile.readline()
            for line in infile:
                row = line.strip().replace('"', '').split(sep)
                row_index = row_map[row[0]]
                for j, elem in enumerate(row[1:]):
                    row_membership[row_index][j] = int(elem)

    return seed


def make_file_column_seeder(filename, sep=' '):
    def seed(matrix, row_membership, num_clusters,
             num_clusters_per_column):
        column_map = {name: idx
                      for idx, name in enumerate(matrix.column_names)}
        column_members = [[] for _ in range(len(matrix.column_names))]

        with open(filename) as infile:
            header = infile.readline()
            for line in infile:
                row = line.strip().replace('"', '').split(sep)
                col_index = column_map[row[0]]
                column_members[col_index] = map(int, row[1:])

        return column_members

    return seed


def make_db_row_seeder(session):
    def seed(row_membership, matrix):
        row_map = {name: idx
                   for idx, name in enumerate(matrix.row_names)}
        iteration = session.query(func.max(cm2db.RowMember.iteration))
        row_clusters = defaultdict(list)
        for row_memb in session.query(cm2db.RowMember).filter(cm2db.RowMember.iteration == iteration):
            row_clusters[row_memb.row_name.name].append(row_memb.cluster)

        # copy memberships
        for row_name in row_clusters.keys():
            #05-07-15 added '.upper' to fix a name mismatch in resume
            if row_name in row_map.keys():
                cur_map = row_map[row_name]
            elif row_name.upper() in row_map.keys():
                cur_map = row_map[row_name.upper()]
            elif row_name.lower() in row_map.keys():
                cur_map = row_map[row_name.lower()]
            else:
                continue

            for i, cluster in enumerate(row_clusters[row_name]):
                #logging.info('row_name = %s, cur_map = %d, cluster = %d, i = %d', row_name, cur_map, cluster, i)
                if i < len(row_membership[cur_map]):
                    row_membership[cur_map][i] = cluster
                else:   #A resumed job that has been finalized may have genes in additional clusters
                    logging.info('Making row_membership for gene %s bigger than %d', cur_map, len(row_membership[cur_map]))
                    row_membership[cur_map].append(cluster)

    return seed


def make_db_column_seeder(session):
    def seed(matrix, row_membership, num_clusters,
             num_clusters_per_column):
        column_map = {name: idx
                      for idx, name in enumerate(matrix.column_names)}

        iteration = session.query(func.max(cm2db.ColumnMember.iteration))
        col_clusters = defaultdict(list)
        for col_memb in session.query(cm2db.ColumnMember).filter(cm2db.ColumnMember.iteration == iteration):
            col_clusters[col_memb.column_name.name].append(col_memb.cluster)

        result = [[0] * num_clusters_per_column for col in matrix.column_names]
        for col_name in matrix.column_names:
            cur_map = column_map[col_name]
            first_expand = True
            for i, cluster in enumerate(col_clusters[col_name]):
                if i < len(result[cur_map]):
                    result[cur_map][i] = cluster
                else:   #A resumed job that has been finalized may have extra columns in clusters
                    if first_expand == True:
                        logging.info('Making col_membership for condition %s bigger, probably due to resuming a finalized run', cur_map)
                        first_expand = False
                    # Make sure that a cluster cannot have the same condition twice
                    if not cluster in result[cur_map]:
                        if 0 in result[cur_map]:
                            result[cur_map][result[cur_map].index(0)] = cluster
                        else:
                            result[cur_map].append(cluster)

        return result
    return seed

def fuzzify(membership, row_scores, column_scores, num_iterations, iteration_result,
            add_fuzz):
    """Provide an iteration-specific fuzzification"""
    if add_fuzz == 'none':
        logging.debug('DO NOT FUZZIFY !!')
        return row_scores, column_scores

    # truth table maps from add_fuzz parameter to where fuzz should be added
    fuzz_vals = {'both': (True, True), 'rows': (True, False), 'columns': (False, True)}
    fuzz_rows, fuzz_cols = fuzz_vals[add_fuzz]

    iteration = iteration_result['iteration']
    #logging.debug("__fuzzify(), setup...")
    #start_time = util.current_millis()
    fuzzy_coeff = old_fuzzy_coefficient(iteration, num_iterations)
    iteration_result['fuzzy-coeff'] = fuzzy_coeff

    if fuzz_rows:
        num_row_fuzzy_values = row_scores.num_rows * row_scores.num_columns
        row_score_values = row_scores.values
        row_sd_values = []

        # iterate the row names directly
        row_names = row_scores.row_names
        for col in xrange(row_scores.num_columns):
            cluster_rows = membership.rows_for_cluster(col + 1)
            for row in xrange(row_scores.num_rows):
                if row_names[row] in cluster_rows:
                    row_sd_values.append(row_score_values[row, col])

        # Note: If there are no non-NaN values in row_sd_values, row_rnorm
        # will have all NaNs
        row_rnorm = util.sd_rnorm(row_sd_values, num_row_fuzzy_values, fuzzy_coeff)
        row_score_values += np.array(row_rnorm).reshape(row_scores.num_rows,
                                                        row_scores.num_columns)

    if fuzz_cols:
        num_col_fuzzy_values = column_scores.num_rows * column_scores.num_columns
        col_score_values = column_scores.values
        col_sd_values = []
        row_names = column_scores.row_names
        for col in xrange(column_scores.num_columns):
            cluster_cols = membership.columns_for_cluster(col + 1)
            for row in xrange(column_scores.num_rows):
                if row_names[row] in cluster_cols:
                    col_sd_values.append(col_score_values[row, col])

        # Note: If there are no non-NaN values in col_sd_values, col_rnorm
        # will have all NaNs
        col_rnorm = util.sd_rnorm(col_sd_values, num_col_fuzzy_values, fuzzy_coeff)
        # add fuzzy values to the row/column scores
        col_score_values += np.array(col_rnorm).reshape(column_scores.num_rows,
                                                        column_scores.num_columns)

    #elapsed = util.current_millis() - start_time
    #logging.debug("fuzzify() finished in %f s.", elapsed / 1000.0)
    return row_scores, column_scores
