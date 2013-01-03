"""iteration_test.py - integration tests for scoring iterations

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import meme
import motif
import unittest
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

KEGG_FILE_PATH = 'testdata/KEGG_taxonomy'
GO_FILE_PATH = 'testdata/proteome2taxid'
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
CACHE_DIR = 'cache'


class IterationTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """This class tests a Halo setup which was extracted from iteration 49
    on an R reference run"""

    def __read_members(self):
        """reads the membership from a fixed set based on halo_ratios5.tsv
        run with the R version"""
        with open('testdata/row_memb-49.tsv') as member_mapfile:
            member_lines = member_mapfile.readlines()
        row_members = {}
        for line in member_lines:
            row = line.strip().split(' ')
            gene = row[0].replace("\"", "")
            for col in range(1, len(row)):
                if gene not in row_members.keys():
                    row_members[gene] = []
                row_members[gene].append(int(row[col]))

        with open('testdata/col_memb-49.tsv') as member_mapfile:
            member_lines = member_mapfile.readlines()
        column_members = {}
        for line in member_lines:
            row = line.strip().split(' ')
            cond = row[0].replace("\"", "")
            for col in range(1, len(row)):
                if cond not in column_members.keys():
                    column_members[cond] = []
                column_members[cond].append(int(row[col]))

        return memb.ClusterMembership(row_members, column_members,
                                      self.config_params)

    def setUp(self):  # pylint; disable-msg=C0103
        """test fixture"""
        self.search_distances = {'upstream': (-20, 150)}
        self.scan_distances = {'upstream': (-30, 250)}

        matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
        infile = util.DelimitedFile.read('halo_ratios5.tsv', has_header=True, quote='\"')
        self.ratio_matrix = matrix_factory.create_from(infile)
        self.organism = make_halo(self.ratio_matrix, self.search_distances, self.scan_distances)
        self.config_params = {'memb.min_cluster_rows_allowed': 3,
                              'memb.max_cluster_rows_allowed': 70,
                              'multiprocessing': False,
                              'memb.clusters_per_row': 2,
                              'memb.clusters_per_col': int(round(43 * 2.0 / 3.0)),
                              'num_clusters': 43,
                              'output_dir': 'out',
                              'num_iterations': 2000}
        self.membership = self.__read_members()  # relies on config_params
        self.iteration_result = { 'iteration': 51 }

    def test_row_scoring(self):
        # tests the row scoring by itself, which combines scoring and fixing
        # extreme values
        row_scoring = microarray.RowScoringFunction(
            self.membership, self.ratio_matrix,
            scaling_func=lambda iteration: 6.0,
            config_params=self.config_params)
        rowscores = row_scoring.compute(self.iteration_result)
        ref_rowscores = read_matrix('testdata/rowscores_fixed.tsv')
        self.assertTrue(check_matrix_values(rowscores, ref_rowscores))

    def test_col_scoring(self):
        # tests the column scoring by itself
        colscoring = scoring.ColumnScoringFunction(
            self.membership, self.ratio_matrix, config_params=self.config_params)
        colscores = colscoring.compute(self.iteration_result)
        ref_colscores = read_matrix('testdata/colscores_fixed.tsv')
        self.assertTrue(check_matrix_values(colscores, ref_colscores))

    def test_net_scoring(self):
        #tests the network scoring by itself#
        network_scaling_fun = scoring.get_default_network_scaling(2000)
        network_scoring = nw.ScoringFunction(self.organism,
                                             self.membership,
                                             self.ratio_matrix,
                                             scaling_func=network_scaling_fun,
                                             run_in_iteration=lambda x: True,
                                             config_params=self.config_params)
        netscores = network_scoring.compute(self.iteration_result).sorted_by_row_name()
        ref_netscores = read_matrix('testdata/netscores_fixed.tsv')
        self.assertTrue(check_matrix_values(netscores, ref_netscores))

    def test_motif_scoring(self):
        meme_suite = meme.MemeSuite430()
        sequence_filters = [
            motif.unique_filter,
            motif.get_remove_low_complexity_filter(meme_suite),
            motif.get_remove_atgs_filter(self.search_distances['upstream'])]
        motif_scaling_fun = scoring.get_default_motif_scaling(2000)
        motif_scoring = motif.MemeScoringFunction(
            self.organism,
            self.membership,
            self.ratio_matrix,
            meme_suite,
            sequence_filters=sequence_filters,
            scaling_func=motif_scaling_fun,
            num_motif_func=motif.default_nmotif_fun,
            update_in_iteration=lambda x: True,
            motif_in_iteration=lambda x: True,
            config_params=self.config_params)
        motscores = motif_scoring.compute(self.iteration_result).sorted_by_row_name()
        ref_motscores = read_matrix('testdata/motscores_fixed.tsv')
        self.assertTrue(check_matrix_values(motscores, ref_motscores))

    def test_scoring_combine(self):
        # tests the scoring combiner
        ref_netscores = read_matrix('testdata/netscores_fixed.tsv')
        ref_motscores = read_matrix('testdata/motscores_fixed.tsv')
        ref_rowscores = read_matrix('testdata/rowscores_fixed.tsv')

        class DummyNetworkScoring(scoring.ScoringFunctionBase):
            def __init__(self):
                scaling_fun = scoring.get_default_network_scaling(2000)
                scoring.ScoringFunctionBase.__init__(self, None, None, scaling_fun)

            def compute(self, iteration_result, ref_matrix=None):
                return ref_netscores

        class DummyMotifScoring(scoring.ScoringFunctionBase):
            def __init__(self):
                scaling_fun = scoring.get_default_motif_scaling(2000)
                scoring.ScoringFunctionBase.__init__(self, None, None, scaling_fun)

            def compute(self, iteration_result, ref_matrix=None):
                return ref_motscores

        class DummyRowScoring(scoring.ScoringFunctionBase):
            def __init__(self):
                scoring.ScoringFunctionBase.__init__(self, None, None,
                                                     lambda iteration: 6.0)

            def compute(self, iteration_result, ref_matrix=None):
                return ref_rowscores

        row_scoring_functions = [DummyRowScoring(), DummyMotifScoring(), DummyNetworkScoring()]
        combiner = scoring.ScoringFunctionCombiner(self.membership,
                                                   row_scoring_functions,
                                                   config_params={'quantile_normalize': True},
                                                   log_subresults=False)
        scores = combiner.compute(self.iteration_result)
        ref_scores = read_matrix('testdata/combined_scores.tsv')
        # note that the rounding error get pretty large here !!!
        self.assertTrue(check_matrix_values(scores, ref_scores))

    def test_quantile_normalize(self):
        row_scores = read_matrix('testdata/rowscores_fixed.tsv')
        mot_scores = read_matrix('testdata/motscores_fixed.tsv')
        net_scores = read_matrix('testdata/netscores_fixed.tsv')

        ref_rowscores = read_matrix('testdata/rowscores_qnorm.tsv')
        ref_motscores = read_matrix('testdata/motscores_qnorm.tsv')
        ref_netscores = read_matrix('testdata/netscores_qnorm.tsv')

        in_matrices = [row_scores, mot_scores, net_scores]
        # scaling for cluster 49
        scalings = [6.0, 0.033355570380253496, 0.016677785190126748]
        result = dm.quantile_normalize_scores(in_matrices, scalings)
        self.assertTrue(check_matrix_values(result[0], ref_rowscores))
        self.assertTrue(check_matrix_values(result[1], ref_motscores))
        self.assertTrue(check_matrix_values(result[2], ref_netscores))

    def test_density_scores(self):
        # tests density score computation
        row_scores = read_matrix('testdata/combined_scores.tsv')
        col_scores = read_matrix('testdata/combined_colscores.tsv')
        ref_rowscores = read_matrix('testdata/density_rowscores.tsv')
        ref_colscores = read_matrix('testdata/density_colscores.tsv')
        rds, cds = memb.get_density_scores(self.membership, row_scores, col_scores)
        self.assertTrue(check_matrix_values(rds, ref_rowscores, eps=1e-11))
        self.assertTrue(check_matrix_values(cds, ref_colscores, eps=1e-11))

    def test_size_compensation(self):
        # tests the size compensation
        row_scores = read_matrix('testdata/density_rowscores.tsv')
        col_scores = read_matrix('testdata/density_colscores.tsv')
        ref_rowscores = read_matrix('testdata/sizecomp_rowscores.tsv')
        ref_colscores = read_matrix('testdata/sizecomp_colscores.tsv')
        memb.compensate_size(self.membership, self.ratio_matrix,
                             row_scores, col_scores)
        self.assertTrue(check_matrix_values(row_scores, ref_rowscores, eps=1e-11))
        self.assertTrue(check_matrix_values(col_scores, ref_colscores, eps=1e-11))


def read_matrix(filename):
    """reads a matrix file"""
    infile = util.DelimitedFile.read(filename, has_header=True, quote='\"')
    return dm.DataMatrixFactory([]).create_from(infile).sorted_by_row_name()


EPS = 1.0e-5

def check_matrix_values(matrix1, matrix2, eps=EPS):
    result = True
    for row in range(matrix1.num_rows):
        for col in range(matrix1.num_columns):
            diff = abs(matrix1.values[row][col] - matrix2.values[row][col])
            #print diff
            if diff > eps:
                print "Value mismatch at (%s, cluster %d): %f != %f (diff = %f)" % (
                    matrix1.row_names[row], col + 1,
                    matrix1.values[row][col],
                    matrix2.values[row][col], diff)
                result = False
    return result


def make_halo(ratio_matrix, search_distances, scan_distances):
    """returns the organism object to work on"""
    keggfile = util.DelimitedFile.read(KEGG_FILE_PATH, comment='#')
    gofile = util.DelimitedFile.read(GO_FILE_PATH)
    rsatdb = rsat.RsatDatabase(RSAT_BASE_URL, CACHE_DIR        )
    mo_db = microbes_online.MicrobesOnline()
    stringfile = 'string_links_64091.tab'

    nw_factories = []
    if stringfile != None:
        nw_factories.append(stringdb.get_network_factory2(stringfile, 0.5))
    else:
        logging.warn("no STRING file specified !")

    nw_factories.append(microbes_online.get_network_factory(
            mo_db, max_operon_size=ratio_matrix.num_rows / 20, weight=0.5))

    org_factory = org.MicrobeFactory(org.make_kegg_code_mapper(keggfile),
                                     org.make_rsat_organism_mapper(rsatdb),
                                     org.make_go_taxonomy_mapper(gofile),
                                     mo_db,
                                     nw_factories)

    return org_factory.create('hal', search_distances, scan_distances)


LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'

if __name__ == '__main__':
    logging.basicConfig(format=LOG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
            IterationTest))
    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
