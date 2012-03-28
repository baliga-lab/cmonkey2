# vi: sw=4 ts=4 et:
"""tps.py - cMonkey tps specific module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import os
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


CHECKPOINT_INTERVAL = 100
CHECKPOINT_FILE = None
CACHE_DIR = 'tpscache'

THESAURUS_FILE = 'tps/tps.synonyms.gz'

NUM_CLUSTERS = 250
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 2000
NETWORK_SCORE_INTERVAL = 7
MOTIF_SCORE_INTERVAL = 10
MAX_CLUSTER_ROWS = 100

SEQUENCE_TYPES = ['upstream']
SEARCH_DISTANCES = {'upstream': (0, 450),'downstream': (0,600)}
SCAN_DISTANCES = {'upstream': (0, 450),'downstream': (0,600)}
UPSTREAM_FILE = 'tps/tps.upstream.-400.50.csv'
DOWNSTREAM_FILE = 'tps/tps.downstream.-100.500.csv'
SEQ_FILENAMES = {'upstream': UPSTREAM_FILE, 'downstream': DOWNSTREAM_FILE }
MAX_MOTIF_WIDTH = 20
STRING_LINKS = 'tps/string_links.tps.tab'

"""these are the default meme iterations ("meme.iters") in the R version"""
MOTIF_ITERS = range( 600, 1200, 100 ) + \
             range( 1250, 1500, 50 ) + \
             range( 1525, 1800, 25 ) + \
             range( 1810, max( NUM_ITERATIONS, 1820 ) + 10 )

mode = 'normal'
#mode = 'debug'
#mode = 'short'
if mode == 'debug':
    NUM_ITERATIONS = 200
    MOTIF_ITERS = [5,100,200]
    NETWORK_SCORE_INTERVAL = 5
    NUM_CLUSTERS = 100

if mode == 'short':
    NUM_ITERATIONS = 500
    MOTIF_ITERS = [100] + range(200,500,50)

def motif_iterations(iteration):
    return iteration in MOTIF_ITERS

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
        self.CHECKPOINT_INTERVAL = CHECKPOINT_INTERVAL

    def organism(self):
        if self.__organism == None:
            self.__organism = self.make_tps()
        return self.__organism

    def make_tps(self):
        """returns a tps organism object"""
        nw_factories = [stringdb.get_network_factory2(STRING_LINKS, 1.0)]
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
        background_file_upstream = meme.global_background_file(
            self.organism(), self.ratio_matrix.row_names(), 'upstream',
            use_revcomp=True)
        background_file_downstream = meme.global_background_file(
            self.organism(), self.ratio_matrix.row_names(), 'downstream',
            use_revcomp=True)
        meme_suite_upstream = meme.MemeSuite430(
            max_width=MAX_MOTIF_WIDTH,
            background_file=background_file_upstream)
        meme_suite_downstream = meme.MemeSuite430(
            max_width=MAX_MOTIF_WIDTH,
            background_file=background_file_downstream)

        upstream_motif_scoring = motif.MemeScoringFunction(
            self.organism(),
            self.membership(),
            self.ratio_matrix,
            meme_suite_upstream,
            seqtype='upstream',
            sequence_filters=sequence_filters,
            pvalue_filter=motif.MinPValueFilter(-20.0),
            weight_func=lambda iteration: 0.0,
            run_in_iteration=motif_iterations,
            config_params=self.config_params)

        downstream_motif_scoring = motif.MemeScoringFunction(
            self.organism(),
            self.membership(),
            self.ratio_matrix,
            meme_suite_downstream,
            seqtype='downstream',
            sequence_filters=sequence_filters,
            pvalue_filter=motif.MinPValueFilter(-20.0),
            weight_func=lambda iteration: 0.0,
            run_in_iteration=motif_iterations,
            config_params=self.config_params)

# NOTE: it looks like Weeder scoring isn't fully implemented yet? At this time, all other run configs seem to be passing MemeSuite430 into the WeederScoringFuncion
#        weeder_scoring = motif.WeederScoringFunction(
#            self.organism(), self.membership(), self.ratio_matrix,
#            meme_suite_downstream, 'downstream',
#            pvalue_filter=motif.MinPValueFilter(-20.0),
#            weight_func=lambda iteration: 0.0,
#            run_in_iteration=motif_iterations,
#            config_params=self.config_params)

        network_scoring = nw.ScoringFunction(self.organism(),
                                             self.membership(),
                                             self.ratio_matrix,
                                             network_iterations,
                                             scoring.default_network_iterations,
                                             config_params=self.config_params)

        motif_combiner = scoring.ScoringFunctionCombiner(
            self.membership(),
            [upstream_motif_scoring, downstream_motif_scoring],
            weight_func=lambda iteration: 0.5)

        return scoring.ScoringFunctionCombiner(
            self.membership(),
            [row_scoring, motif_combiner, network_scoring])


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
        ratios = sys.argv[1]
        if mode == 'debug': ratios = 'tps_ratios.1000x4.tsv'
        infile = util.DelimitedFile.read(ratios, has_header=True, quote='\"')
        matrix = matrix_factory.create_from(infile)
        cmonkey_run = TpsCMonkeyRun('tps', matrix, NUM_CLUSTERS)
        if os.path.exists(CHECKPOINT_FILE): cmonkey_run.run_from_checkpoint(CHECKPOINT_FILE)
        else: cmonkey_run.run()
