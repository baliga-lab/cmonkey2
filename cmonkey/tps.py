# vi: sw=4 ts=4 et:
"""tps.py - cMonkey tps specific module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import organism
import scoring
import microarray
import stringdb
import network as nw
import motif
import meme
import datamatrix as dm
import membership as memb
import util
import cmonkey_run


CACHE_DIR = 'tpscache'

THESAURUS_FILE = 'tps/tps.synonyms.gz'

NUM_CLUSTERS = 250
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 2000
NETWORK_SCORE_INTERVAL = 7
MOTIF_SCORE_INTERVAL = 10
MAX_CLUSTER_ROWS = 100

SEQUENCE_TYPES = ['upstream']
SEARCH_DISTANCES = {'upstream': (0, 400)}
SCAN_DISTANCES = {'upstream': (0, 400)}
PROM_SEQFILE = 'tps/tps.upstream.-350.50.csv'
SEQ_FILENAMES = {'upstream': PROM_SEQFILE}
MAX_MOTIF_WIDTH = 20
STRING_LINKS = 'tps/string_links.tps.tab'

"""these are the default meme iterations ("meme.iters") in the R version"""
MEME_ITERS = range( 600, 1200, 100 ) + \
             range( 1250, 1500, 50 ) + \
             range( 1525, 1800, 25 ) + \
             range( 1810, max( NUM_ITERATIONS, 1820 ) + 10 )

mode = 'normal'
#mode = 'debug'
#mode = 'short'
if mode == 'debug':
    NUM_ITERATIONS = 200
    MEME_ITERS = [2,100,200]
    NETWORK_SCORE_INTERVAL = 5
if mode == 'short':
    NUM_ITERATIONS = 500
    MEME_ITERS = [100] + range(200,500,50)

def meme_iterations(iteration):
    return iteration in MEME_ITERS

def network_iterations(iteration):
    return iteration > 0 and iteration % NETWORK_SCORE_INTERVAL == 0


class TpsCMonkeyRun(cmonkey_run.CMonkeyRun):

    def __init__(self, organism_code, ratio_matrix, num_clusters):
        cmonkey_run.CMonkeyRun.__init__(self, organism_code, ratio_matrix, num_clusters)
        self.__organism = None
        self['cache_dir'] = CACHE_DIR
        self['sequence_types'] = SEQUENCE_TYPES
        self['search_distances'] = SEARCH_DISTANCES
        self['scan_distances'] = SCAN_DISTANCES

    def organism(self):
        if self.__organism == None:
            self.__organism = self.make_tps()
        return self.__organism

    def make_tps(self):
        """returns a tps organism object"""
        nw_factories = [stringdb.get_network_factory2(STRING_LINKS)]
        return organism.GenericOrganism('tps', THESAURUS_FILE, nw_factories,
                                        seq_filenames=SEQ_FILENAMES,
                                        search_distances=SEARCH_DISTANCES,
                                        scan_distances=SCAN_DISTANCES)

    def make_row_scoring(self):
        """returns the row scoring function"""
        row_scoring = microarray.RowScoringFunction(
            self.membership(), self.ratio_matrix,
            lambda iteration: ROW_WEIGHT,
            config_params=self.config_params)

        sequence_filters = []
        background_file_prom = meme.global_background_file(
            self.organism(), self.ratio_matrix.row_names(), 'upstream',
            use_revcomp=True)
        meme_suite_prom = meme.MemeSuite430(
            max_width=MAX_MOTIF_WIDTH,
            background_file=background_file_prom)

        motif_scoring = motif.MemeScoringFunction(
            self.organism(),
            self.membership(),
            self.ratio_matrix,
            meme_suite_prom,
            seqtype='upstream',
            sequence_filters=sequence_filters,
            pvalue_filter=motif.MinPValueFilter(-20.0),
            weight_func=lambda iteration: 0.0,
            run_in_iteration=meme_iterations,
            config_params=self.config_params)

        network_scoring = nw.ScoringFunction(self.organism(),
                                             self.membership(),
                                             self.ratio_matrix,
                                             network_iterations,
                                             scoring.default_network_iterations,
                                             config_params=self.config_params)

        return scoring.ScoringFunctionCombiner(
            self.membership(),
            [row_scoring, motif_scoring, network_scoring])


if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) < 2:
        print('Usage: ./rembrandt.sh <ratio-file> [checkpoint-file]')
    else:
        if len(sys.argv) > 2:
            CHECKPOINT_FILE = sys.argv[2]
        matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
        infile = util.DelimitedFile.read(sys.argv[1], has_header=True, quote='\"')
        matrix = matrix_factory.create_from(infile)
        cmonkey_run = TpsCMonkeyRun('tps', matrix, NUM_CLUSTERS)
        cmonkey_run.run()
