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


LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
CMONKEY_VERSION = '4.0'
KEGG_FILE_PATH = 'testdata/KEGG_taxonomy'
GO_FILE_PATH = 'testdata/proteome2taxid'
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
CACHE_DIR = 'cache'
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 2000
NETWORK_SCORE_INTERVAL = 7
MOTIF_SCORE_INTERVAL = 10

# used to select sequences and MEME
SEARCH_DISTANCES = {'upstream': (-20, 150)}
# used for background distribution and MAST
SCAN_DISTANCES = {'upstream': (-30, 250)}


def run_cmonkey():
    """init of the cMonkey system"""
    logging.basicConfig(format=LOG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    if not os.path.exists(CACHE_DIR):
        os.mkdir(CACHE_DIR)

    matrix = read_matrix(sys.argv[1]).sorted_by_row_name()
    logging.info("Normalized input matrix has %d rows and %d columns:",
                 matrix.num_rows(), matrix.num_columns())

    organism = make_organism(sys.argv[2], matrix)
    membership = make_membership(matrix)
    gene_scoring = make_gene_scoring_func(organism, membership, matrix)
    cond_scoring = microarray.ColumnScoringFunction(membership, matrix)

    for iteration in range(NUM_ITERATIONS):
        iterate(membership, matrix, gene_scoring, cond_scoring,
                iteration, NUM_ITERATIONS)
    print "Done !!!!"
    print "cluster\t# rows"
    for cluster in range(1, membership.num_clusters() + 1):
        print "%d\t%d" % (cluster, membership.num_row_members(cluster))
    #print membership


def make_gene_scoring_func(organism, membership, matrix):
    """setup the gene-related scoring functions here
    each object in this array supports the method
    compute(organism, membership, matrix) and returns
    a DataMatrix(genes x cluster)
    """
    row_scoring = microarray.RowScoringFunction(membership, matrix,
                                                lambda iteration: ROW_WEIGHT)

    meme_suite = meme.MemeSuite430()
    sequence_filters = [
        motif.unique_filter,
        motif.get_remove_low_complexity_filter(meme_suite),
        motif.get_remove_atgs_filter(SEARCH_DISTANCES['upstream'])]

    motif_scoring = motif.ScoringFunction(
        organism,
        membership,
        matrix,
        meme_suite,
        sequence_filters=sequence_filters,
        pvalue_filter=motif.make_min_value_filter(-20.0),
        weight_func=lambda iteration: 0.0,
        interval=MOTIF_SCORE_INTERVAL)

    network_scoring = nw.ScoringFunction(organism, membership, matrix,
                                         lambda iteration: 0.0,
                                         NETWORK_SCORE_INTERVAL)

    return memb.ScoringFunctionCombiner([row_scoring, motif_scoring,
                                         network_scoring])


def read_matrix(filename):
    """reads the data matrix from a file"""
    matrix_factory = dm.DataMatrixFactory(
        [dm.nochange_filter, dm.center_scale_filter])
    infile = util.DelimitedFile.read(filename, has_header=True)
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
        stringdb.get_network_factory2(stringdb.STRING_FILE2),
        microbes_online.get_network_factory(
            mo_db, max_operon_size=matrix.num_rows() / 20)]
    org_factory = org.MicrobeFactory(org.make_kegg_code_mapper(keggfile),
                                     org.make_rsat_organism_mapper(rsatdb),
                                     org.make_go_taxonomy_mapper(gofile),
                                     mo_db,
                                     nw_factories)
    return org_factory.create(organism_code, SEARCH_DISTANCES, SCAN_DISTANCES)


def make_membership(matrix):
    """returns a seeded membership"""
    # We are using a fake row seed here in order to have reproducible,
    # deterministic results for development. Since the column seed is
    # deterministic and dependent on the row seed, we can use the real
    # column seed implementation here.
    # We currently assume the halo_ratios5.tsv file as input
    fake_row_membership_seed = util.DelimitedFileMapper(
        util.DelimitedFile.read('clusters.tsv', has_header=False), 0, 1)
    return memb.ClusterMembership.create(
        matrix.sorted_by_row_name(),
        fake_seed_row_memberships(fake_row_membership_seed),
        microarray.seed_column_members)


def iterate(membership, matrix, gene_scoring_func, cond_scoring_func,
            iteration, num_iterations):
    """one iteration of the algorithm"""
    logging.info("Iteration # %d", iteration)
    membership.update(matrix,
                      gene_scoring_func.compute(iteration),
                      cond_scoring_func.compute(iteration),
                      iteration, num_iterations)

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
