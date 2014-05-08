"""postproc_test.py - integration tests for post iteration adjustments

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import meme
import motif
import unittest
import xmlrunner
import util
import rsat
import organism as org
import microbes_online
import stringdb
import datamatrix as dm
import membership as memb
import microarray
import scoring
import network as nw
import logging
import sys
import testutil


class PostprocTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """This class tests a Halo setup which was extracted from iteration 49
    on an R reference run"""

    def __read_members(self):
        """reads the membership from a fixed set based on halo_ratios5.tsv
        run with the R version"""
        with open('testdata/row_membs43-before-postproc.tsv') as member_mapfile:
            member_lines = member_mapfile.readlines()
        row_members = {}
        for line in member_lines:
            row = line.strip().split('\t')
            gene = row[0].replace("\"", "")
            for col in range(1, len(row)):
                if gene not in row_members.keys():
                    row_members[gene] = []
                row_members[gene].append(int(row[col]))

        with open('testdata/col_membs43-before-postproc.tsv') as member_mapfile:
            member_lines = member_mapfile.readlines()
        column_members = {}
        for line in member_lines:
            row = line.strip().split('\t')
            cond = row[0].replace("\"", "")
            for col in range(1, len(row)):
                if cond not in column_members.keys():
                    column_members[cond] = []
                column_members[cond].append(int(row[col]))

        return memb.OrigMembership(sorted(row_members.keys()),
                                   sorted(column_members.keys()),
                                   row_members, column_members,
                                   self.config_params)

    def setUp(self):  # pylint; disable-msg=C0103
        """test fixture"""
        self.search_distances = {'upstream': (-20, 150)}
        self.scan_distances = {'upstream': (-30, 250)}

        matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
        infile = util.read_dfile('example_data/hal/halo_ratios5.tsv',
                                 has_header=True, quote='\"')
        self.ratio_matrix = matrix_factory.create_from(infile)
        self.organism = testutil.make_halo(self.search_distances, self.scan_distances,
                                           self.ratio_matrix)
        self.config_params = {'memb.min_cluster_rows_allowed': 3,
                              'memb.max_cluster_rows_allowed': 70,
                              'multiprocessing': False,
                              'memb.clusters_per_row': 2,
                              'memb.clusters_per_col': int(round(43 * 2.0 / 3.0)),
                              'num_clusters': 43,
                              'num_iterations': 2000}
        self.membership = self.__read_members()  # relies on config_params
        self.iteration_result = { 'iteration': 51 }

    def test_adjust_cluster(self):
        """tests the row scoring by itself, which combines scoring and fixing
        extreme values"""
        rowscores = read_matrix('testdata/rowscores-43-before-postproc.tsv')
        assign1 = memb.adjust_cluster(self.membership,
                                      1, rowscores, cutoff=0.33,
                                      limit=100)
        self.assertEquals(0, len(assign1))
        assign2 = memb.adjust_cluster(self.membership,
                                      2, rowscores, cutoff=0.33,
                                      limit=100)
        self.assertEquals(1, len(assign2))
        self.assertEquals(2, assign2['VNG6210G'])
        assign5 = memb.adjust_cluster(self.membership,
                                      5, rowscores, cutoff=0.33,
                                      limit=100)
        self.assertEquals(3, len(assign5))
        self.assertEquals(5, assign5['VNG1182H'])
        self.assertEquals(5, assign5['VNG2259C'])
        self.assertEquals(5, assign5['VNG0719G'])

    def test_post_adjust(self):
        """tests the row scoring by itself, which combines scoring and fixing
        extreme values"""
        rowscores = read_matrix('testdata/rowscores-43-before-postproc.tsv')
        memb.postadjust(self.membership, rowscores)


def read_matrix(filename):
    """reads a matrix file"""
    infile = util.read_dfile(filename, has_header=True, quote='\"')
    return dm.DataMatrixFactory([]).create_from(infile).sorted_by_row_name()


EPS = 1.0e-5

def check_matrix_values(matrix1, matrix2, eps=EPS):
    result = True
    for row in range(matrix1.num_rows()):
        for col in range(matrix1.num_columns()):
            diff = abs(matrix1[row][col] - matrix2[row][col])
            #print diff
            if diff > eps:
                print "Value mismatch at (%s, cluster %d): %f != %f (diff = %f)" % (
                    matrix1.row_name(row), col + 1, matrix1[row][col], matrix2[row][col], diff)
                result = False
    return result


LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'

if __name__ == '__main__':
    logging.basicConfig(format=LOG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
            PostprocTest))
    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(SUITE))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
