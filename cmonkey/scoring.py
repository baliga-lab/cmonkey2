"""scoring.py - cMonkey scoring base classes

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'

import microarray
import logging
import os
import datamatrix as dm


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

    def weight(self, iteration):
        """returns the weight for the specified iteration"""
        return self.__weight_func(iteration)


class ScoringFunctionCombiner:
    """Taking advantage of the composite pattern, this combiner function
    exposes the basic interface of a scoring function in order to
    allow for nested scoring functions as they are used in the motif
    scoring
    """
    def __init__(self, scoring_functions, weight_func=None):
        """creates a combiner instance"""
        self.__scoring_functions = scoring_functions

    def compute(self, iteration):
        """compute scores for one iteration"""
        result_matrices = []
        score_weights = []
        for scoring_function in self.__scoring_functions:
            matrix = scoring_function.compute(iteration)
            if matrix != None:
                result_matrices.append(matrix)
                score_weights.append(scoring_function.weight(iteration))

        if len(result_matrices) > 0:
            result_matrices = dm.quantile_normalize_scores(result_matrices,
                                                           score_weights)
        combined_score = (result_matrices[0] *
                          self.__scoring_functions[0].weight(iteration))
        for index in range(1, len(result_matrices)):
            combined_score += (
                result_matrices[index] *
                self.__scoring_functions[index].weight(iteration))
        return combined_score

    def weight(self, iteration):
        """returns the weight for the specified iteration"""
        return self.__weight_func(iteration)


class ConfigurationBase:
    """configuration base class"""

    def __init__(self, organism_code, matrix_filename,
                 num_iterations, cache_dir):
        """create instance"""
        logging.basicConfig(format=LOG_FORMAT,
                            datefmt='%Y-%m-%d %H:%M:%S',
                            level=logging.DEBUG)
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)
        self.__cache_dir = cache_dir
        self.__matrix_filename = matrix_filename
        self.__organism_code = organism_code
        self.__num_iterations = num_iterations

        self.__matrix = None
        self.__membership = None
        self.__organism = None
        self.__row_scoring = None
        self.__column_scoring = None

    def organism_code(self):
        """returns the organism code"""
        return self.__organism_code

    def num_iterations(self):
        """returns the number of iterations"""
        return self.__num_iterations

    def cache_dir(self):
        """returns the cache directory"""
        return self.__cache_dir

    def matrix(self):
        """returns the input matrix"""
        if self.__matrix == None:
            self.__matrix = self.read_matrix(
                self.__matrix_filename).sorted_by_row_name()
        return self.__matrix

    def membership(self):
        """returns the seeded membership"""
        if self.__membership == None:
            self.__membership = self.make_membership()
        return self.__membership

    def organism(self):
        """returns the organism object to work on"""
        if self.__organism == None:
            self.__organism = self.make_organism()
        return self.__organism

    def row_scoring(self):
        """returns the row scoring function"""
        if self.__row_scoring == None:
            self.__row_scoring = self.make_row_scoring()
        return self.__row_scoring

    def column_scoring(self):
        """returns the column scoring function"""
        if self.__column_scoring == None:
            self.__column_scoring = microarray.ColumnScoringFunction(
                self.membership(), self.matrix())
        return self.__column_scoring

    def make_membership(self):
        """implement in derived class"""
        pass

    def read_matrix(self, filename):
        """implement in derived class"""
        pass

    def make_organism(self):
        """implement in derived class"""
        pass

    def make_row_scoring(self):
        """implement in derived class"""
        pass
