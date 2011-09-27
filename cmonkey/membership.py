"""membership.py - cMonkey cluster membership functionality
This module captures the microarray-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import datamatrix as dm
import math
import util

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
    def __init__(self, num_clusters, row_is_member_of, column_is_member_of):
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
        self.__row_is_member_of = row_is_member_of
        self.__column_is_member_of = column_is_member_of
        self.__cluster_row_members = create_cluster_to_names_map(
            self.__row_is_member_of)
        self.__cluster_column_members = create_cluster_to_names_map(
            self.__column_is_member_of)

    # pylint: disable-msg=R0913
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
        return ClusterMembership(num_clusters, row_is_member_of,
                                 col_is_member_of)

    def num_clusters(self):
        """returns the number of clusters"""
        return self.__num_clusters

    def clusters_for_row(self, row_name):
        """determine the clusters for the specified row"""
        return self.__row_is_member_of[row_name]

    def clusters_for_column(self, column_name):
        """determine the clusters for the specified column"""
        return self.__column_is_member_of[column_name]

    def rows_for_cluster(self, cluster):
        """determine the rows that belong to a cluster"""
        return sorted(self.__cluster_row_members[cluster])

    def columns_for_cluster(self, cluster):
        """determine the rows that belong to a cluster"""
        return sorted(self.__cluster_column_members[cluster])

    def is_row_member_of(self, row_name, cluster):
        """determines whether a certain row is member of a cluster"""
        return row_name in self.rows_for_cluster(cluster)


class Membership:
    """Algorithms for cluster membership"""
    def __init__(self):
        """this class has only class methods"""
        pass

    @classmethod
    def map_to_is_member_matrix(cls, membership_matrix, kcluster):
        """maps a matrix containing row/column numbers to a true/false
        matrix by checking all values i in the range [1, kcluster] for
        containment in each row of the membership matrix.
        Example: mat =  [1 2] kcluster: 3
                        [2 3]
                 result = [t f] => 1 is in     [1 2], but not in [2 3]
                          [t t] => 2 is in     [1 2], and        [2 3]
                          [f t] => 3 is not in [1 2], but in     [2 3]
        """
        result = []
        for i in range(1, kcluster + 1):
            result_row = []
            result.append(result_row)
            for matrix_row in membership_matrix:
                result_row.append(i in matrix_row)

        return result


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


def get_density_scores(membership, row_scores, col_scores, num_clusters):
    """We can't really implement density scores at the moment,
    there seems to be no equivalent to R's density() and approx()
    in scipy"""
    rscore_range = abs(row_scores.max() - row_scores.min())
    rowscore_bandwidth = max(rscore_range / 100.0, 0.001)
    rd_scores = dm.DataMatrix(row_scores.num_rows(),
                              row_scores.num_columns(),
                              row_scores.row_names(),
                              row_scores.column_names())
    for cluster in range(1, num_clusters + 1):
        rr_scores = get_rr_scores(membership, row_scores, rowscore_bandwidth,
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
        cc_scores = get_cc_scores(membership, col_scores, colscore_bandwidth,
                                  cluster)
        for row in range(col_scores.num_rows()):
            cd_scores[row][cluster - 1] = cc_scores[row]
    return (rd_scores, cd_scores)


def get_rr_scores(membership, rowscores, bandwidth, cluster):
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


def get_cc_scores(membership, scores, bandwidth, cluster):
    """calculate the density scores for the given column score values in the
    specified cluster"""
    cluster_columns = membership.columns_for_cluster(cluster)

    row_names = scores.row_names()
    score_indexes = [row_names.index(name) for name in cluster_columns]
    kscores = scores.column_values(cluster - 1)
    cluster_scores = [kscores[index] for index in score_indexes]
    return util.density(kscores, cluster_scores, bandwidth,
                        min(kscores) - 1, max(kscores) + 1)
