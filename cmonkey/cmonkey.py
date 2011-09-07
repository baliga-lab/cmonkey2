"""cmonkey.py - cMonkey top-level module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import util
import datamatrix as dm
import microarray
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


class CMonkey:  # pylint: disable-msg=R0902
    """
    The cMonkey object controls the overall execution of the cMonkey
    algorithm.
    This top-level class takes configuration and inputs to provide
    them for the actual execution
    organism: Organism instance to use for computation
    ratio_matrices: a list of DataMatrix instances containing gene expressions
    configuration is a dictionary with the following keys:
    'num_iterations'   - # of iterations
    'clusters_per_row' - # of clusters/row
    ...
    """

    def __init__(self, organism, ratio_matrices, config=None):
        """create a cMonkey object
        ratio_matrices: a MatrixCollection object containing gene expression
                        values
        configuration: a dictionary of configuration values
        """
        self.__finished = False
        self.__organism = organism
        self.__ratio_matrices = ratio_matrices
        self.__configuration = self.__init_configuration(config)
        self.__networks = None

        # we might move these to a result class
        self.__row_scores = None
        self.__col_scores = None
        self.__gene_weights = None

    def __init_configuration(self, config):
        """sets meaningful defaults in the configuration"""
        configuration = config or {}
        configuration.setdefault('num_iterations',       2)  # 2000 for now
        configuration.setdefault('biclusters_per_gene',  2)
        configuration.setdefault('num_biclusters',
                                 self.__compute_num_biclusters(configuration))
        configuration.setdefault('row.scaling',          6)
        configuration.setdefault('seed_method.rows',     'kmeans')
        configuration.setdefault('seed_method.cols',     'best')
        return configuration

    def __compute_num_biclusters(self, config):
        """computes the number of biclusters to optimize"""
        return int(round(self.__ratio_matrices.num_unique_rows() *
                         float(config['biclusters_per_gene']) /
                         AVG_CLUSTER_SIZE))

    def __num_biclusters(self):
        """returns the number of biclusters"""
        return self.__configuration['num_biclusters']

    def finished(self):
        """returns true if the run was finished"""
        return self.__finished

    def run(self):
        """start a run"""
        self.__retrieve_networks()
        logging.info("# Networks read: %d", len(self.__networks))

        # self.__seed_clusters()
        current_iteration = 0

        logging.info("# iterations: %d", self.__num_iterations())

        while current_iteration < self.__num_iterations():
            self.__iterate()
            current_iteration += 1
        self.__finished = True

    def __input_gene_names(self):
        """returns the unique gene names used in the input matrices"""
        return self.__ratio_matrices[0].row_names()

    def __retrieve_networks(self):
        """retrieves the networks provided by the organism object and
        possibly other sources, doing some normalization if necessary"""
        self.__networks = self.__organism.networks()
        max_score = 0
        for network in self.__networks:
            logging.info("Network with %d edges", network.num_edges())
            nw_total = network.total_score()
            if nw_total > max_score:
                max_score = nw_total
        for network in self.__networks:
            network.normalize_scores_to(max_score)

    def __num_iterations(self):
        """configured number of iterations"""
        return self.__configuration.get('num_iterations', 2000)

    def __iterate(self):
        """iteration step in cMonkey
        This is run over and over until no improvements can be achieved"""
        pass


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
    # 43 clusters, so it follows:
    # n.clust.per.row = 2
    # n.clust.per.col = k.clust * 2/3 => 43 * 2/3 => 29
    membership = microarray.ClusterMembership.create(
        matrix.sorted_by_row_name(), 43, 2,
        29, fake_seed_row_memberships(fake_row_membership_seed),
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
    #microarray.compute_row_scores(membership, matrix, NUM_CLUSTERS)
    #cscores = microarray.compute_column_scores(membership, matrix,
    #                                           NUM_CLUSTERS)
    #print cscores

    # 2. compute motif scores
    #meme_suite = meme.MemeSuite430()
    #distance = motif.DISTANCE_UPSTREAM_SEARCH
    #motif.compute_scores(meme_suite, organism, membership, used_seqs,
    #                     distance,
    #                     [motif.unique_filter,
    #                      motif.get_remove_low_complexity_filter(meme_suite),
    #                      motif.remove_atgs_filter])

    # running the algorithm in the CMonkey object is obsolete
    #algorithm = CMonkey(organism, dm.DataMatrixCollection([matrix]))
    #algorithm.run()
    # 3. compute network scores
    #print "# USED SEQS: %d" % len(matrix.row_names())
    #cluster_genes = ['VNG1691G', 'VNG1699C', 'VNG1709G', 'VNG1713G',
    #                 'VNG1715G', 'VNG1718G']
    networks = retrieve_networks(organism)
    for network in networks:
        logging.info("Compute scores for network '%s'", network.name())
        for cluster_num in range(1, NUM_CLUSTERS + 1):
            cluster_genes = sorted(membership.rows_for_cluster(cluster_num))
            network_scores = nw.compute_network_scores(network,
                                                       cluster_genes,
                                                       matrix.row_names())
            for gene, score in network_scores:
                print "%s[%d] - %s\t%f" % (network.name(), cluster_num, gene, score)


def retrieve_networks(organism):
    """retrieves the networks provided by the organism object and
    possibly other sources, doing some normalization if necessary"""
    networks = organism.networks()
    max_score = 0
    for network in networks:
        logging.info("Network '%s' with %d edges", network.name(),
                     network.num_edges())
        nw_total = network.total_score()
        if nw_total > max_score:
            max_score = nw_total
    for network in networks:
        network.normalize_scores_to(max_score)
    return networks


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
