"""cmonkey.py - cMonkey top-level module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import logging
import scoring
import microarray
import datamatrix as dm
import membership as memb
import util

CMONKEY_VERSION = '4.0'
CHECKPOINT_INTERVAL = 3
CHECKPOINT_FILE = None

NUM_ITERATIONS = 2000
CACHE_DIR = 'leishmania_cache'
NUM_CLUSTERS = 267
ROW_WEIGHT = 6.0

THESAURUS_FILE = 'leishmania_data/synonymThesaurus.csv.gz'

SEQUENCE_TYPES = ['upstream', 'p3utr']
SEARCH_DISTANCES = {'upstream': (0, 700), 'p3utr': (0, 831)}
SCAN_DISTANCES = {'upstream': (0, 700), 'p3utr': (0, 831)}
PROM_SEQFILE = 'leishmania_data/promoterSeqs_Mus_musculus_RA.csv'
P3UTR_SEQFILE = 'leishmania_data/p3utrSeqs_Mus_musculus_RA.csv.gz'
SEQ_FILENAMES = {'upstream': PROM_SEQFILE, 'p3utr': P3UTR_SEQFILE}


class CMonkeyConfiguration(scoring.ConfigurationBase):
    """Leishmania-specific configuration class"""
    def __init__(self, config_params,
                 checkpoint_file=None):
        """create instance"""
        scoring.ConfigurationBase.__init__(self, config_params,
                                           checkpoint_file)

    @classmethod
    def create(cls, matrix_filename, checkpoint_file=None):
        """creates an initialized instance"""
        params = (scoring.ConfigurationBuilder().
                  with_organism('mmu').
                  with_matrix_filenames([matrix_filename]).
                  with_num_iterations(NUM_ITERATIONS).
                  with_cache_dir(CACHE_DIR).
                  with_num_clusters(NUM_CLUSTERS).
                  with_sequence_types(SEQUENCE_TYPES).
                  with_search_distances(SEARCH_DISTANCES).
                  with_scan_distances(SCAN_DISTANCES).
                  build())
        return cls(params, checkpoint_file)

    def read_matrix(self, filename):
        """returns the matrix"""
        #return read_matrix(filename)
        matrix_factory = dm.DataMatrixFactory(
            [dm.nochange_filter, dm.center_scale_filter])
        infile = util.DelimitedFile.read(filename, has_header=True,
                                         quote='\"')
        return matrix_factory.create_from(infile)

    def make_membership(self):
        """returns the seeded membership"""
        num_clusters = self.config_params[memb.KEY_NUM_CLUSTERS]
        return memb.ClusterMembership.create(
            self.matrix().sorted_by_row_name(),
            memb.make_kmeans_row_seeder(num_clusters),
            microarray.seed_column_members,
            self.config_params)

    def make_organism(self):
        """returns a configured organism object"""
        nw_factories = [stringdb.get_network_factory2('leishmania_data/mouse_links_preprocessed.csv', sep=';')]
        return organism.GenericOrganism('mmu', THESAURUS_FILE, nw_factories,
                                        seq_filenames=SEQ_FILENAMES,
                                        search_distances=SEARCH_DISTANCES,
                                        scan_distances=SCAN_DISTANCES)


    def make_row_scoring(self):
        """returns the row scoring function"""
        row_scoring = microarray.RowScoringFunction(
            self.membership(), self.matrix(),
            lambda iteration: ROW_WEIGHT,
            config_params=self.config_params)

        return scoring.ScoringFunctionCombiner(
            self.membership(), [row_scoring])


def run_cmonkey(config):
    """init of the cMonkey system"""
    gene_scoring = config.row_scoring()
    cond_scoring = config.column_scoring()

    for iteration in range(config.start_iteration(),
                           config.num_iterations()):
        logging.info("Iteration # %d", iteration)
        config.membership().update(config.matrix(),
                                   gene_scoring.compute(iteration),
                                   cond_scoring.compute(iteration),
                                   iteration, config.num_iterations())
        if iteration > 0 and  iteration % CHECKPOINT_INTERVAL == 0:
            config.save_checkpoint_data(iteration)
    print "Done !!!!"
    print "cluster\t# rows"
    for cluster in range(1, config.membership().num_clusters() + 1):
        print "%d\t%d" % (cluster,
                          config.membership().num_row_members(cluster))
    #print membership


if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) < 2:
        print('Usage: ./leishmania.sh <ratio-file> [checkpoint-file]')
    else:
        if len(sys.argv) > 2:
            CHECKPOINT_FILE = sys.argv[2]
        
        conf = CMonkeyConfiguration.create(
            sys.argv[1], checkpoint_file=CHECKPOINT_FILE)
        run_cmonkey(conf)
