"""microbial_setup.py - cMonkey module for microbe-specific configuration

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import util
import logging
import datamatrix as dm
import microarray
import scoring
import organism as org
import meme
import motif
import network as nw
import microbes_online
import stringdb
import rsat
import membership as memb

KEGG_FILE_PATH = 'testdata/KEGG_taxonomy'
GO_FILE_PATH = 'testdata/proteome2taxid'
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
CACHE_DIR = 'cache'
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 2000
NETWORK_SCORE_INTERVAL = 7
MOTIF_SCORE_INTERVAL = 10
NUM_CLUSTERS = 43

SEQUENCE_TYPES = ['upstream']
# used to select sequences and MEME
SEARCH_DISTANCES = {'upstream': (-20, 150)}
# used for background distribution and MAST
SCAN_DISTANCES = {'upstream': (-30, 250)}


class CMonkeyConfiguration(scoring.ConfigurationBase):
    """Microbe-specific configuration class"""

    def __init__(self, config_params, checkpoint_file=None):
        """create instance"""
        scoring.ConfigurationBase.__init__(self, config_params,
                                           checkpoint_file)

    @classmethod
    def create(cls, organism_code, matrix_filename,
               checkpoint_file=None):
        """Creates an initialized instance"""
        params = (scoring.ConfigurationBuilder().
                  with_organism(organism_code).
                  with_num_iterations(NUM_ITERATIONS).
                  with_matrix_filenames([matrix_filename]).
                  with_cache_dir(CACHE_DIR).
                  with_num_clusters(NUM_CLUSTERS).
                  with_sequence_types(SEQUENCE_TYPES).
                  with_search_distances(SEARCH_DISTANCES).
                  with_scan_distances(SCAN_DISTANCES).
                  build())
        print "CONFIG: ", params
        return cls(params, checkpoint_file)

    def read_matrix(self, filename):
        """returns the matrix"""
        return read_matrix(filename)

    def make_membership(self):
        """returns the seeded membership"""
        return make_membership(self.matrix())

    def make_organism(self):
        """returns the organism object to work on"""
        return make_organism(self.organism_code(), self.matrix())

    def row_scoring(self):
        """returns the row scoring function"""
        return make_gene_scoring_func(self.organism(),
                                      self.membership(),
                                      self.matrix())


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

    motif_scoring = motif.MemeScoringFunction(
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

    return scoring.ScoringFunctionCombiner([row_scoring, motif_scoring,
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
