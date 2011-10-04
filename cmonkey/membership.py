"""membership.py - cMonkey cluster membership functionality
This module captures the microarray-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import datamatrix as dm
import math
import util
import random
import logging
import sys


# Default values for membership creation
NUM_CLUSTERS = 43
NUM_CLUSTERS_PER_ROW = 2
NUM_CLUSTERS_PER_COL = int(round(NUM_CLUSTERS * 2.0 / 3.0))
MIN_CLUSTER_ROWS_ALLOWED = 3
MAX_CLUSTER_ROWS_ALLOWED = 70

PROBABILITY_SEEING_ROW_CHANGE = 0.5
PROBABILITY_SEEING_COLUMN_CHANGE = 1.0
MAX_CHANGES_PER_ROW = 1
MAX_CHANGES_PER_COLUMN = 5


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
    def __init__(self,
                 row_is_member_of,
                 column_is_member_of,
                 num_clusters,
                 num_clusters_per_row,
                 num_clusters_per_col,
                 probability_seeing_row_change,
                 probability_seeing_col_change,
                 max_changes_per_row,
                 max_changes_per_column):
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

        self.__num_clusters = num_clusters
        self.__num_clusters_per_row = num_clusters_per_row
        self.__num_clusters_per_col = num_clusters_per_col
        self.__probability_seeing_row_change = probability_seeing_row_change
        self.__probability_seeing_column_change = probability_seeing_col_change
        self.__max_changes_per_row = max_changes_per_row
        self.__max_changes_per_column = max_changes_per_column

        self.__row_is_member_of = row_is_member_of
        self.__column_is_member_of = column_is_member_of
        self.__cluster_row_members = create_cluster_to_names_map(
            self.__row_is_member_of)
        self.__cluster_column_members = create_cluster_to_names_map(
            self.__column_is_member_of)

    # pylint: disable-msg=R0913
    @classmethod
    def create(cls, data_matrix,
               seed_row_memberships,
               seed_column_memberships,
               num_clusters=NUM_CLUSTERS,
               num_clusters_per_row=NUM_CLUSTERS_PER_ROW,
               num_clusters_per_column=NUM_CLUSTERS_PER_COL,
               probability_seeing_row_change=PROBABILITY_SEEING_ROW_CHANGE,
               probability_seeing_col_change=PROBABILITY_SEEING_COLUMN_CHANGE,
               max_changes_per_row=MAX_CHANGES_PER_ROW,
               max_changes_per_col=MAX_CHANGES_PER_COLUMN):
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
        return ClusterMembership(row_is_member_of,
                                 col_is_member_of,
                                 num_clusters,
                                 num_clusters_per_row,
                                 num_clusters_per_column,
                                 probability_seeing_row_change,
                                 probability_seeing_col_change,
                                 max_changes_per_row,
                                 max_changes_per_col)

    def num_clusters(self):
        """returns the number of clusters"""
        return self.__num_clusters

    def num_clusters_per_row(self):
        """returns the number of clusters per row"""
        return self.__num_clusters_per_row

    def num_clusters_per_column(self):
        """returns the number of clusters per row"""
        return self.__num_clusters_per_col

    def clusters_for_row(self, row_name):
        """determine the clusters for the specified row"""
        if row_name in self.__row_is_member_of:
            return self.__row_is_member_of[row_name]
        else:
            return []

    def num_clusters_for_row(self, row_name):
        """returns the number of clusters for the row"""
        return len(self.clusters_for_row(row_name))

    def clusters_for_column(self, column_name):
        """determine the clusters for the specified column"""
        if column_name in self.__column_is_member_of:
            return self.__column_is_member_of[column_name]
        else:
            return []

    def num_clusters_for_column(self, column_name):
        """returns the number of clusters for the column"""
        return len(self.clusters_for_column(column_name))

    def rows_for_cluster(self, cluster):
        """determine the rows that belong to a cluster"""
        return sorted(self.__cluster_row_members[cluster])

    def columns_for_cluster(self, cluster):
        """determine the rows that belong to a cluster"""
        return sorted(self.__cluster_column_members[cluster])

    def is_row_member_of(self, row_name, cluster):
        """determines whether a certain row is member of a cluster"""
        return row_name in self.rows_for_cluster(cluster)

    def is_column_member_of(self, column_name, cluster):
        """determines whether a certain column is member of a cluster"""
        return column_name in self.columns_for_cluster(cluster)

    def num_row_members(self, cluster):
        """returns the number of row members in the specified cluster"""
        if cluster in self.__cluster_row_members:
            return len(self.__cluster_row_members[cluster])
        else:
            return 0

    def num_column_members(self, cluster):
        """returns the number of row members in the specified cluster"""
        if cluster in self.__cluster_column_members:
            return len(self.__cluster_column_members[cluster])
        else:
            return 0

    def is_row_in_clusters(self, row, clusters):
        """returns true if the specified row is in all spefied clusters"""
        for cluster in clusters:
            if not self.is_row_member_of(row, cluster):
                return False
        return True

    def is_column_in_clusters(self, col, clusters):
        """returns true if the specified row is in all spefied clusters"""
        for cluster in clusters:
            if not self.is_column_member_of(col, cluster):
                return False
        return True

    def add_cluster_to_row(self, row, cluster):
        """checked adding of a row to a cluster"""
        if self.num_clusters_for_row(row) >= self.__num_clusters_per_row:
            raise Exception(("add_row_to_cluster() - exceeded clusters/row " +
                            "limit for row: '%s'" % str(row)))
        self.__add_cluster_to_row(row, cluster)

    def __add_cluster_to_row(self, row, cluster):
        """adds the specified row as a member to the cluster. Unchecked
        version, without checking limits"""
        if not row in self.__row_is_member_of:
            self.__row_is_member_of[row] = []
        if not cluster in self.__cluster_row_members:
            self.__cluster_row_members[cluster] = []

        clusters = self.__row_is_member_of[row]
        rows = self.__cluster_row_members[cluster]

        if cluster not in clusters:
            logging.info("ROW %s -> CLUSTER %d", row, cluster)
            clusters.append(cluster)
        else:
            pass
            #logging.warn("cluster %s already associated with %s",
            #             str(cluster), str(row))
        if row not in rows:
            rows.append(row)

    def remove_cluster_from_row(self, row, cluster):
        """removes a cluster from the list of associated clusters for a row"""
        if row in self.__row_is_member_of:
            clusters = self.__row_is_member_of[row]
            clusters.remove(cluster)

        if cluster in self.__cluster_row_members:
            rows = self.__cluster_row_members[cluster]
            rows.remove(row)

    def replace_row_cluster(self, row, cluster, replacement):
        """replaces a cluster in the list of clusters for a row"""
        self.remove_cluster_from_row(row, cluster)
        self.add_cluster_to_row(row, replacement)

    def add_cluster_to_column(self, column, cluster):
        """checked adding of a column to a cluster"""
        if self.num_clusters_for_column(column) >= self.__num_clusters_per_col:
            raise Exception(("add_col_to_cluster() - exceeded clusters/col " +
                            "limit for col: '%s'" % str(column)))
        self.__add_cluster_to_column(column, cluster)

    def __add_cluster_to_column(self, column, cluster):
        """adds the specified column as a member to the cluster. Unchecked
        version, without checking limits"""
        if not column in self.__column_is_member_of:
            self.__column_is_member_of[column] = []
        if not cluster in self.__cluster_column_members:
            self.__cluster_column_members[cluster] = []

        clusters = self.__column_is_member_of[column]
        columns = self.__cluster_column_members[cluster]

        if cluster not in clusters:
            logging.info("COL %s -> CLUSTER %d", column, cluster)
            clusters.append(cluster)
        else:
            #logging.warn("cluster %s already associated with %s",
            #             str(cluster), str(column))
            pass
        if columns not in columns:
            columns.append(column)

    def remove_cluster_from_column(self, column, cluster):
        """removes a cluster from the list of associated clusters for a column"""
        if cluster in self.__cluster_column_members:
            columns = self.__cluster_column_members[cluster]
            columns.remove(column)

        if column in self.__column_is_member_of:
            clusters = self.__column_is_member_of[column]
            clusters.remove(cluster)
 
    def replace_column_cluster(self, column, cluster, replacement):
        """replaces a cluster in the list of clusters for a column"""
        self.remove_cluster_from_column(column, cluster)
        self.add_cluster_to_column(column, replacement)

    def __repr__(self):
        """returns the string representation of memberships"""
        result = "ROW MEMBERS:\n"
        result += repr(self.__cluster_row_members)
        result += "\n\nCOLUMN MEMBERS:\n"
        result += repr(self.__cluster_column_members)
        return result

    def update(self, matrix, row_scores, column_scores):
        """top-level update method"""
        # TODO: Fuzzify scores (can't be reproduced 1:1 to the R version)
        rd_scores, cd_scores = _get_density_scores(self, row_scores,
                                                   column_scores)
        _compensate_size(self, matrix, rd_scores, cd_scores)
        self._update_memberships(rd_scores, cd_scores)

    def _update_memberships(self, rd_scores, cd_scores):
        """update memberships according to rd_scores and cd_scores"""

        def get_best_clusters(scores, num_per_cluster):
            """retrieve the best scored gene clusters from the given
            row/column score matrix"""
            result = {}
            for row in range(scores.num_rows()):
                row_values = scores.row_values(row)
                ranked_scores = sorted(row_values, reverse=True)
                rowname = scores.row_names()[row]
                result[rowname] = []
                for index in range(num_per_cluster):
                    result[rowname].append(row_values.index(
                            ranked_scores[index]) + 1)
            return result

        def seeing_change(prob):
            """returns true if the update is seeing the change"""
            return prob >= 1.0 or random.uniform(0.0, 1.0) > prob

        def add_cluster_to_row(row, cluster):
            """ Ways to add a member to a cluster:
            1. if the number of members is less than the allowed, simply add
            2. if there is a conflict, replace a gene with a lower score in the
               scores matrix
            """
            if self.num_clusters_for_row(row) < self.__num_clusters_per_row:
                self.add_cluster_to_row(row, cluster)
            else:
                replace_lowest_scoring_row_member(row, cluster)

        def replace_lowest_scoring_row_member(row, cluster):
            current_clusters = self.__row_is_member_of[row]
            min_score = sys.maxint
            min_cluster = None
            member_names = rd_scores.row_names()
            member_index = member_names.index(row)

            for current_cluster in current_clusters:
                if rd_scores[member_index][current_cluster - 1] < min_score:
                    min_score = rd_scores[member_index][current_cluster - 1]
                    min_cluster = current_cluster
            self.remove_cluster_from_row(row, min_cluster)
            self.add_cluster_to_row(row, cluster)

        def add_cluster_to_col(col, cluster):
            """adds a column to a cluster"""
            if self.num_clusters_for_column(col) < self.__num_clusters_per_col:
                self.add_cluster_to_column(col, cluster)
            else:
                replace_lowest_scoring_col_member(col, cluster)

        def replace_lowest_scoring_col_member(column, cluster):
            current_clusters = self.__column_is_member_of[column]
            min_score = sys.maxint
            min_cluster = None
            member_names = cd_scores.row_names()
            member_index = member_names.index(column)

            for current_cluster in current_clusters:
                if cd_scores[member_index][current_cluster - 1] < min_score:
                    min_score = cd_scores[member_index][current_cluster - 1]
                    min_cluster = current_cluster
            self.remove_cluster_from_column(column, min_cluster)
            self.add_cluster_to_column(column, cluster)

        def update_for(scores,
                       num_per_cluster,
                       probability_seeing_change,
                       is_in_all_clusters, max_changes,
                       add_member_to_cluster):
            """generically updating row/column memberships according to
            rd_scores/cd_scores"""
            best_clusters = get_best_clusters(scores, num_per_cluster)
            for row in range(scores.num_rows()):
                rowname = scores.row_names()[row]
                best_members = best_clusters[rowname]
                if (not is_in_all_clusters(rowname, best_members) and
                    seeing_change(probability_seeing_change)):
                    for change in range(max_changes):
                        add_member_to_cluster(rowname, best_members[change])

        update_for(rd_scores,
                   self.__num_clusters_per_row,
                   self.__probability_seeing_row_change,
                   self.is_row_in_clusters, self.__max_changes_per_row,
                   add_cluster_to_row)
        update_for(cd_scores,
                   self.__num_clusters_per_col,
                   self.__probability_seeing_column_change,
                   self.is_column_in_clusters, self.__max_changes_per_column,
                   add_cluster_to_col)


class ScoringFunctionBase:
    """Base class for scoring functions"""

    def __init__(self, membership, matrix, weight_func):
        """creates a function instance"""
        self.__membership = membership
        self.__matrix = matrix
        self.__weight_func = weight_func

    def membership(self):
        """returns this function's membership object"""
        return self.__membership

    def matrix(self):
        """returns this function's matrix object"""
        return self.__matrix

    def compute(self, iteration):
        """general compute method, iteration is the 0-based iteration number"""
        raise Exception("please implement me")

    def num_clusters(self):
        """returns the number of clusters"""
        return self.__membership.num_clusters()

    def gene_names(self):
        """returns the gene names"""
        return self.__matrix.row_names()

    def rows_for_cluster(self, cluster):
        """returns the rows for the specified cluster"""
        return self.__membership.rows_for_cluster(cluster)

    def apply_weight(self, result, iteration):
        """applies the stored weight function to the result
        Note: we might be able to incorporate this into the result object"""
        if (self.__weight_func != None):
            return result.multiply_by(self.__weight_func(iteration))


def _get_density_scores(membership, row_scores, col_scores):
    """We can't really implement density scores at the moment,
    there seems to be no equivalent to R's density() and approx()
    in scipy"""
    num_clusters = membership.num_clusters()
    rscore_range = abs(row_scores.max() - row_scores.min())
    rowscore_bandwidth = max(rscore_range / 100.0, 0.001)
    rd_scores = dm.DataMatrix(row_scores.num_rows(),
                              row_scores.num_columns(),
                              row_scores.row_names(),
                              row_scores.column_names())
    for cluster in range(1, num_clusters + 1):
        rr_scores = _get_rr_scores(membership, row_scores, rowscore_bandwidth,
                                   cluster)
        for row in range(row_scores.num_rows()):
            rd_scores[row][cluster - 1] = rr_scores[row]

    cscore_range = abs(col_scores.max() - col_scores.min())
    colscore_bandwidth = max(cscore_range / 100.0, 0.001)
    cd_scores = dm.DataMatrix(col_scores.num_rows(),
                              col_scores.num_columns(),
                              col_scores.row_names(),
                              col_scores.column_names())
    for cluster in range(1, num_clusters + 1):
        cc_scores = _get_cc_scores(membership, col_scores, colscore_bandwidth,
                                   cluster)
        for row in range(col_scores.num_rows()):
            cd_scores[row][cluster - 1] = cc_scores[row]
    return (rd_scores, cd_scores)


def _get_rr_scores(membership, rowscores, bandwidth, cluster):
    """calculate the density scores for the given row score values in the
    specified cluster"""
    def bwscale(value):
        """standard bandwidth scaling function for row scores"""
        return math.exp(-value / 10.0) * 10.0

    cluster_rows = membership.rows_for_cluster(cluster)

    row_names = rowscores.row_names()
    score_indexes = [row_names.index(name) for name in cluster_rows]
    kscores = rowscores.column_values(cluster - 1)
    cluster_scores = [kscores[index] for index in score_indexes]
    cluster_bandwidth = bandwidth * bwscale(len(cluster_rows))

    return util.density(kscores, cluster_scores, cluster_bandwidth,
                        min(kscores) - 1, max(kscores) + 1)


def _get_cc_scores(membership, scores, bandwidth, cluster):
    """calculate the density scores for the given column score values in the
    specified cluster"""
    cluster_columns = membership.columns_for_cluster(cluster)

    row_names = scores.row_names()
    score_indexes = [row_names.index(name) for name in cluster_columns]
    kscores = scores.column_values(cluster - 1)
    cluster_scores = [kscores[index] for index in score_indexes]
    return util.density(kscores, cluster_scores, bandwidth,
                        min(kscores) - 1, max(kscores) + 1)


def _compensate_size(membership, matrix, rd_scores, cd_scores):
    """size compensation function"""
    def compensate_dim_size(size, dimsize, clusters_per_dim, num_clusters):
        """compensate size for a dimension"""
        return math.exp(-size / (dimsize * clusters_per_dim) / num_clusters)

    def compensate_row_size(size):
        """compensation function for row dimension"""
        return compensate_dim_size(size,
                                   matrix.num_rows(),
                                   membership.num_clusters_per_row(),
                                   membership.num_clusters())

    def compensate_column_size(size):
        """compensation function for column dimension"""
        return compensate_dim_size(size,
                                   matrix.num_columns(),
                                   membership.num_clusters_per_column(),
                                   membership.num_clusters())

    def compensate_rows(cluster):
        """compensate density scores for row dimension"""
        num_rowmembers = membership.num_row_members(cluster)
        if num_rowmembers > 0:
            rd_scores.multiply_column_by(
                cluster - 1, compensate_row_size(num_rowmembers))
        else:
            rd_scores.multiply_column_by(
                cluster - 1, compensate_row_size(MIN_CLUSTER_ROWS_ALLOWED))

    def compensate_columns(cluster):
        """compensate density scores for column dimension"""
        num_colmembers = membership.num_column_members(cluster)
        if num_colmembers > 0:
            cd_scores.multiply_column_by(
                cluster - 1, compensate_column_size(num_colmembers))
        else:
            cd_scores.multiply_column_by(
                cluster - 1,
                compensate_column_size(matrix.num_columns() / 10.0))

    num_clusters = membership.num_clusters()
    for cluster in range(1, num_clusters + 1):
        compensate_rows(cluster)
        compensate_columns(cluster)
