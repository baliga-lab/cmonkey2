"""microarray_test.py - unit test module for microarray module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import datamatrix as dm
import util
import membership as memb
import microarray as ma
import numpy


class ComputeArrayScoresTest(unittest.TestCase):
    """compute_row_scores"""
    def __read_members(self):
        with open('testdata/row_membership.tsv') as member_mapfile:
            member_lines = member_mapfile.readlines()
        row_members = {}
        for line in member_lines:
            row = line.strip().split('\t')
            row_members[row[0]] = [int(row[1])]

        with open('testdata/column_membership.tsv') as member_mapfile:
            member_lines = member_mapfile.readlines()
        column_members = {}
        for line in member_lines:
            row = line.strip().split('\t')
            column_members[row[0]] = [int(cluster)
                                      for cluster in row[1].split(':')]
        return memb.ClusterMembership(row_members, column_members,
                                      {'memb.num_clusters': 43})

    def __read_ratios(self):
        dfile = util.DelimitedFile.read('testdata/row_scores_testratios.tsv',
                                        has_header=True)
        return dm.DataMatrixFactory([]).create_from(dfile)

    def __read_rowscores_refresult(self):
        dfile = util.DelimitedFile.read('testdata/row_scores_refresult.tsv',
                                        has_header=True, quote='"')
        return dm.DataMatrixFactory([]).create_from(dfile)

    def __read_colscores_refresult(self):
        dfile = util.DelimitedFile.read('testdata/column_scores_refresult.tsv',
                                        has_header=True, quote='"')
        return dm.DataMatrixFactory([]).create_from(dfile)

    def test_compute_row_scores_multiprocessing(self):
        membership = self.__read_members()
        ratios = self.__read_ratios()
        print "(reading reference row scores...)"
        refresult = self.__read_rowscores_refresult()
        print "(compute my own row scores...)"
        result = ma.compute_row_scores(membership, ratios, 43, True)
        print "(comparing computed with reference results...)"
        self.__compare_with_refresult(refresult, result)

    def test_compute_row_scores_single(self):
        membership = self.__read_members()
        ratios = self.__read_ratios()
        print "(reading reference row scores...)"
        refresult = self.__read_rowscores_refresult()
        print "(compute my own row scores...)"
        result = ma.compute_row_scores(membership, ratios, 43, False)
        print "(comparing computed with reference results...)"
        self.__compare_with_refresult(refresult, result)

    def test_compute_column_scores(self):
        membership = self.__read_members()
        ratios = self.__read_ratios()
        refresult = self.__read_colscores_refresult()
        result = ma.compute_column_scores(membership, ratios, 43)
        self.__compare_with_refresult(refresult, result)

    def __compare_with_refresult(self, refresult, result):
        self.assertEquals(refresult.num_rows(), result.num_rows())
        self.assertEquals(refresult.num_columns(), result.num_columns())
        self.assertTrue((result.row_names() == refresult.row_names()).all())
        for row_index in range(result.num_rows()):
            for col_index in range(result.num_columns()):
                # note that we reduced the comparison's number of places
                # That's because the input matrix was created with some
                # rounding, so we have a slightly higher rounding difference
                self.assertAlmostEquals(refresult[row_index][col_index],
                                        result[row_index][col_index], 3)

    def test_subtract_and_square(self):
        matrix = dm.DataMatrix(2, 3, values=[[11.0, 12.0, 13.0],
                                             [14.0, 15.0, 16.0]])
        result = ma.subtract_and_square(matrix, [5.0, 6.0, 3.0])
        self.assertTrue(([[36.0, 36.0, 100.0],
                          [81.0, 81.0, 169.0]] == result).all())
