"""cmonkey.py - cMonkey top-level module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import util
import datamatrix as dm
import microarray
import membership as memb
import organism as org
import meme
import motif
import network as nw
import microbes_online
import stringdb
import rsat
import sys
import os
import logging
import math

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
CMONKEY_VERSION = '4.0'
AVG_CLUSTER_SIZE = 20
KEGG_FILE_PATH = 'testdata/KEGG_taxonomy'
GO_FILE_PATH = 'testdata/proteome2taxid'
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
CACHE_DIR = 'cache'
NUM_CLUSTERS = 43
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 1

# num_clusters clusters, so it follows for num_clusters = 43:
# num_clusters_per_row = 2
# num_clusters_per_col = num_clusters * 2/3 => 43 * 2/3 => 29
NUM_CLUSTERS_PER_ROW = 2
NUM_CLUSTERS_PER_COL = int(round(NUM_CLUSTERS * 2.0 / 3.0))
MIN_CLUSTER_ROWS_ALLOWED = 3
MAX_CLUSTER_ROWS_ALLOWED = 70


def run_cmonkey():
    """init of the cMonkey system"""
    logging.basicConfig(format=LOG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    if not os.path.exists(CACHE_DIR):
        os.mkdir(CACHE_DIR)

    matrix = read_matrix(sys.argv[1])
    logging.info("Normalized input matrix has %d rows and %d columns:",
                 matrix.num_rows(), matrix.num_columns())

    organism = make_organism(sys.argv[2], matrix)
    membership = make_membership(matrix, NUM_CLUSTERS,
                                 NUM_CLUSTERS_PER_ROW,
                                 NUM_CLUSTERS_PER_COL)

    # microarray scoring
    row_scoring = microarray.RowScoringFunction(membership, matrix,
                                                lambda iteration: ROW_WEIGHT)

    meme_suite = meme.MemeSuite430()
    sequence_filters = [motif.unique_filter,
                        motif.get_remove_low_complexity_filter(meme_suite),
                        motif.remove_atgs_filter]

    motif_scoring = motif.ScoringFunction(organism,
                                          membership,
                                          matrix,
                                          meme_suite,
                                          sequence_filters,
                                          motif.make_min_value_filter(-20.0))

    network_scoring = nw.ScoringFunction(organism, membership, matrix)

    # setup all scoring functions in this array so they are executed
    # one after another.
    # each object in this array supports the method
    # compute(organism, membership, matrix) and returns
    # a DataMatrix(genes x cluster)
    #scoring_funcs = [row_scoring, motif_scoring, network_scoring]
    scoring_funcs = [row_scoring]
    cscoring = microarray.ColumnScoringFunction(membership, matrix)
    for iteration in range(NUM_ITERATIONS):
        iterate(membership, matrix, scoring_funcs, cscoring, iteration)
    print "Done !!!!"


def read_matrix(filename):
    """reads the data matrix from a file"""
    matrix_factory = dm.DataMatrixFactory(
        [dm.nochange_filter, dm.center_scale_filter])
    infile = util.DelimitedFile.read(sys.argv[1], has_header=True)
    # for now, we set a fixed set of clusters assignments
    # that matches halo_ratios5.tsv
    # This is to take out the random component out of the seeding
    # and make it easier to compare results with R-cMonkey
    return matrix_factory.create_from(infile)


def make_organism(organism_code, matrix):
    """create the organism based on the organism code"""
    keggfile = util.DelimitedFile.read(KEGG_FILE_PATH, comment='#')
    gofile = util.DelimitedFile.read(GO_FILE_PATH)
    rsatdb = rsat.RsatDatabase(RSAT_BASE_URL, CACHE_DIR)
    mo_db = microbes_online.MicrobesOnline()
    # note that for the moment, the STRING factory is hardwired to
    # a preprocessed Halobacterium SP file
    nw_factories = [
        microbes_online.get_network_factory(
            mo_db, max_operon_size=matrix.num_rows() / 20),
        stringdb.get_network_factory(stringdb.STRING_FILE2)]
    org_factory = org.OrganismFactory(org.make_kegg_code_mapper(keggfile),
                                      org.make_rsat_organism_mapper(rsatdb),
                                      org.make_go_taxonomy_mapper(gofile),
                                      mo_db,
                                      nw_factories)
    return org_factory.create(organism_code)


def make_membership(matrix, num_clusters, num_clusters_per_row,
                    num_clusters_per_col):
    """returns a seeded membership"""
    logging.info("# CLUSTERS = %d", num_clusters)
    logging.info("# CLUSTERS / ROW = %d", num_clusters_per_row)
    logging.info("# CLUSTERS / COL = %d", num_clusters_per_col)

    # We are using a fake row seed here in order to have reproducible,
    # deterministic results for development. Since the column seed is
    # deterministic and dependent on the row seed, we can use the real
    # column seed implementation here.
    # We currently assume the halo_ratios5.tsv file as input
    fake_row_membership_seed = util.DelimitedFileMapper(
        util.DelimitedFile.read('clusters.tsv', has_header=False), 0, 1)
    return memb.ClusterMembership.create(
        matrix.sorted_by_row_name(),
        num_clusters,
        num_clusters_per_row,
        num_clusters_per_col,
        fake_seed_row_memberships(fake_row_membership_seed),
        microarray.seed_column_members)


def iterate(membership, matrix, scoring_funcs, column_scoring_func, iteration):
    """one iteration of the algorithm"""
    logging.info("Iteration # %d", iteration)
    result_matrices = []
    for score_func in scoring_funcs:
        result_matrices.append(score_func.compute(iteration))
    cscores = column_scoring_func.compute(iteration)

    # TODO: combine (log filter + weight)
    for index in range(len(scoring_funcs)):
        result_matrices[index] = scoring_funcs[index].apply_weight(
            result_matrices[index], iteration)

    # TODO: Fuzzify scores (can't be reproduced 1:1 to the R version)
    # Get density score
    rd_scores, cd_scores = memb.get_density_scores(membership,
                                                   result_matrices[0],
                                                   cscores)
    compensate_size(membership, matrix, rd_scores, cd_scores)


def compensate_size(membership, matrix, rd_scores, cd_scores):
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
                cluster - 1, compensate_column_size(matrix.num_columns() / 10.0))

    num_clusters = membership.num_clusters()
    for cluster in range(1, num_clusters + 1):
        compensate_rows(cluster)
        compensate_columns(cluster)


############################################################
#### Replace with real seeding when everything works
############################################################


def fake_seed_row_memberships(fake_mapper):
    """This method sets the memberships according to a seed that was
    created by running the original cMonkey on halo_ratios5.tsv with
    kmeans row seeding. To compromise on its NP complete behavior,
    kmeans does not always return the same clusters. We bake all random
    components of cMonkey for development to make it possible to compare
    results"""
    def compute(row_membership, _):
        """pseudo-seed with fixed numbers"""
        logging.debug("fake_seed_row_memberships")
        index = 0
        for key in sorted(fake_mapper.keys()):
            row_membership[index][0] = int(fake_mapper[key])
            index += 1
    return compute

if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) <= 2:
        print('Usage: ./run_cmonkey.sh <ratio-file> <organism-code>')
    else:
        run_cmonkey()


__all__ = ['CMonkey', 'Membership']
