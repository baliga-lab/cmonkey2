"""microarray_test.py - unit test module for microarray module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import cmonkey.datamatrix as dm
import cmonkey.util as util
import cmonkey.membership as memb
import cmonkey.microarray as ma
import cmonkey.scoring as scoring
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
        return memb.OrigMembership(sorted(row_members.keys()),
                                   sorted(column_members.keys()),
                                   row_members, column_members,
                                   {'memb.num_clusters': 43,
                                    'memb.clusters_per_row': 2,
                                    'memb.clusters_per_col': 29 })

    def __read_ratios(self):
        dfile = util.read_dfile('testdata/row_scores_testratios.tsv',
                                has_header=True)
        return dm.DataMatrixFactory([]).create_from(dfile, case_sensitive=True)

    def __read_rowscores_refresult(self):
        dfile = util.read_dfile('testdata/row_scores_refresult.tsv',
                                has_header=True, quote='"')
        return dm.DataMatrixFactory([]).create_from(dfile, case_sensitive=True)

    def __read_colscores_refresult(self):
        dfile = util.read_dfile('testdata/column_scores_refresult.tsv',
                                has_header=True, quote='"')
        return dm.DataMatrixFactory([]).create_from(dfile, case_sensitive=True)

    def test_compute_row_scores_multiprocessing(self):
        membership = self.__read_members()
        ratios = self.__read_ratios()
        print("(reading reference row scores...)")
        refresult = self.__read_rowscores_refresult()
        print("(compute my own row scores...)")
        result = ma.compute_row_scores(membership, ratios, 43,
                                       {'multiprocessing': True, 'num_cores': None})
        result.fix_extreme_values()
        print("(comparing computed with reference results...)")
        self.__compare_with_refresult(refresult, result)

    def test_compute_row_scores_single(self):
        membership = self.__read_members()
        ratios = self.__read_ratios()
        print("(reading reference row scores...)")
        refresult = self.__read_rowscores_refresult()
        print("(compute my own row scores...)")
        result = ma.compute_row_scores(membership, ratios, 43,
                                       {'multiprocessing': True, 'num_cores': None})
        result.fix_extreme_values()
        print("(comparing computed with reference results...)")
        self.__compare_with_refresult(refresult, result)

    def test_compute_column_scores(self):
        membership = self.__read_members()
        ratios = self.__read_ratios()
        refresult = self.__read_colscores_refresult()
        result = scoring.compute_column_scores(membership, ratios, 43,
                                               {'multiprocessing': True, 'num_cores': None})
        self.__compare_with_refresult(refresult, result)

    def __compare_with_refresult(self, refresult, result):
        self.assertEquals(refresult.num_rows, result.num_rows)
        self.assertEquals(refresult.num_columns, result.num_columns)
        self.assertEquals(result.row_names, refresult.row_names)
        for row_index in range(result.num_rows):
            for col_index in range(result.num_columns):
                # note that we reduced the comparison's number of places
                # That's because the input matrix was created with some
                # rounding, so we have a slightly higher rounding difference
                self.assertAlmostEquals(refresult.values[row_index][col_index],
                                        result.values[row_index][col_index], 3)
