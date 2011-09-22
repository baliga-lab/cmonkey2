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
AVG_CLUSTER_SIZE = 20
KEGG_FILE_PATH = 'testdata/KEGG_taxonomy'
GO_FILE_PATH = 'testdata/proteome2taxid'
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
CACHE_DIR = 'cache'
NUM_CLUSTERS = 43


def run_cmonkey():
    """init of the cMonkey system"""
    logging.basicConfig(format=LOG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    if not os.path.exists(CACHE_DIR):
        os.mkdir(CACHE_DIR)

    matrix_factory = dm.DataMatrixFactory(
        [dm.nochange_filter, dm.center_scale_filter])
    infile = util.DelimitedFile.read(sys.argv[1], has_header=True)
    # for now, we set a fixed set of clusters assignments
    # that matches halo_ratios5.tsv
    # This is to take out the random component out of the seeding
    # and make it easier to compare results with R-cMonkey
    fake_row_membership_seed = util.DelimitedFileMapper(
        util.DelimitedFile.read('clusters.tsv', has_header=False), 0, 1)
    matrix = matrix_factory.create_from(infile)
    logging.info("Normalized input matrix has %d rows and %d columns:",
                 matrix.num_rows(), matrix.num_columns())
    #print(matrix.sorted_by_row_name())

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
    # We are using a fake row seed here in order to have reproducible,
    # deterministic results for development. Since the column seed is
    # deterministic and dependent on the row seed, we can use the real
    # column seed implementation here.
    # We currently assume the halo_ratios5.tsv file as input, and
    # NUM_CLUSTERS clusters, so it follows:
    # n.clust.per.row = 2
    # n.clust.per.col = k.clust * 2/3 => 43 * 2/3 => 29
    num_clusters_per_row = 2
    num_clusters_per_col = int(round(NUM_CLUSTERS * 2.0 / 3.0))

    logging.info("# CLUSTERS = %d", NUM_CLUSTERS)
    logging.info("# CLUSTERS / ROW = %d", num_clusters_per_row)
    logging.info("# CLUSTERS / COL = %d", num_clusters_per_col)

    membership = memb.ClusterMembership.create(
        matrix.sorted_by_row_name(),
        NUM_CLUSTERS,
        num_clusters_per_row,
        num_clusters_per_col,
        fake_seed_row_memberships(fake_row_membership_seed),
        microarray.seed_column_members)

    organism = org_factory.create(sys.argv[2])

    # precompute the sequences for all genes that are referenced in the
    # input ratios, they are used as a basis to compute the background
    # distribution for every cluster
    used_seqs = organism.sequences_for_genes(sorted(matrix.row_names()),
                                             motif.DISTANCE_UPSTREAM_SCAN,
                                             upstream=True)

    # One iteration
    # 1. compute microarray scores
    # setup all scoring functions in this array so they are executed
    # one after another.
    # each object in this array supports the method
    # compute(organism, membership, matrix) and returns
    # a DataMatrix(genes x cluster)
    row_scoring = microarray.RowScoringFunction(NUM_CLUSTERS)

    #rscores = scoring_algos[0].compute(organism, membership, matrix)
    #rscores = rscores.multiply_by(6.0) # TODO: don't hardcode

    #cscores = microarray.ColumnScoringFunction(NUM_CLUSTERS).compute(
    #    organism, membership, matrix)
    #print cscores

    # 2. compute motif scores
    meme_suite = meme.MemeSuite430()
    sequence_filters = [motif.unique_filter,
                        motif.get_remove_low_complexity_filter(meme_suite),
                        motif.remove_atgs_filter]

    motif_scoring = motif.ScoringFunction(organism,
                                          meme_suite,
                                          motif.DISTANCE_UPSTREAM_SEARCH,
                                          NUM_CLUSTERS,
                                          used_seqs,
                                          sequence_filters,
                                          motif.make_min_value_filter(-20.0))

    # 3. compute network scores
    network_scoring = nw.ScoringFunction(organism, NUM_CLUSTERS)

    scoring_algos = [row_scoring, motif_scoring, network_scoring]
    result_matrices = []
    for score_func in scoring_algos:
        result_matrices.append(score_func.compute(membership, matrix))

    print "Done !!!!"

    # running the algorithm in the CMonkey object is obsolete
    #algorithm = CMonkey(organism, dm.DataMatrixCollection([matrix]))
    #algorithm.run()

    # TODO: Fuzzify scores (can't be reproduced 1:1 to the R version)
    # TODO: Get density score
    # TODO: size compensation


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
