# vi: sw=4 ts=4 et:
"""microbe_custom_run.py - cMonkey specific module for non-standard microbe runs

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import os
import scoring
import microarray
import network as nw
import motif
import meme
import datamatrix as dm
import util
import cmonkey_run


CHECKPOINT_INTERVAL = 100
CHECKPOINT_FILE = None

NUM_CLUSTERS = 250
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 2000
NETWORK_SCORE_INTERVAL = 7
MOTIF_SCORE_INTERVAL = 10
MAX_CLUSTER_ROWS = 100

#SEQUENCE_TYPES = ['upstream']
#SEARCH_DISTANCES = {'upstream': (0, 450)}
#SCAN_DISTANCES = {'upstream': (0, 450)}
#MAX_MOTIF_WIDTH = 20

"""these are the default meme iterations ("meme.iters") in the R version"""
MOTIF_ITERS = range( 600, 1200, 100 ) + \
              range( 1250, 1500, 50 ) + \
              range( 1525, 1800, 25 ) + \
              range( 1810, max( NUM_ITERATIONS, 1820 ), 10 )

mode = 'normal'
mode = 'debug'
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


class MicrobeCustomRun(cmonkey_run.CMonkeyRun):

    def __init__(self, organism_code, ratio_matrix, num_clusters):
        cmonkey_run.CMonkeyRun.__init__(self, organism_code, ratio_matrix, num_clusters)
        self.CHECKPOINT_INTERVAL = CHECKPOINT_INTERVAL

    def make_row_scoring(self):
        """makes a row scoring function on demand"""
        # Default row scoring functions
        row_scoring = microarray.RowScoringFunction(
            self.membership(), self.ratio_matrix,
            scaling_func=lambda iteration: self['row_scaling'],
            config_params=self.config_params)

        meme_suite = meme.MemeSuite481(remove_tempfiles=True)
        sequence_filters = [
            motif.unique_filter,
            motif.get_remove_low_complexity_filter(meme_suite),
            motif.get_remove_atgs_filter(self['search_distances']['upstream'])]

        motif_scoring = motif.MemeScoringFunction(
            self.organism(),
            self.membership(),
            self.ratio_matrix,
            meme_suite,
            sequence_filters=sequence_filters,
            pvalue_filter=motif.MinPValueFilter(-20.0),
            scaling_func=lambda iteration: 1.0,  # TODO
            run_in_iteration=motif_iterations,
            config_params=self.config_params)

        network_scoring = nw.ScoringFunction(self.organism(),
                                             self.membership(),
                                             self.ratio_matrix,
                                             scaling_func=lambda iteration: 0.0,
                                             run_in_iteration=network_iterations,
                                             config_params=self.config_params)

        row_scoring_functions = [row_scoring, motif_scoring, network_scoring]
        return scoring.ScoringFunctionCombiner(self.membership(),
                                               row_scoring_functions,
                                               log_subresults=True)


if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) < 4:
        print('Usage: ./run_cmonkey.sh <organism-code> <ratio-file> <string-file> [checkpoint-file]')
    else:
        if len(sys.argv) > 4: CHECKPOINT_FILE = sys.argv[4]
        matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
        ratios = sys.argv[2]
        infile = util.DelimitedFile.read(ratios, has_header=True, quote='\"')
        matrix = matrix_factory.create_from(infile)
        cmonkey_run = MicrobeCustomRun(sys.argv[1], matrix, NUM_CLUSTERS)
        cmonkey_run['string_file'] = sys.argv[3]
        if CHECKPOINT_FILE and os.path.exists(CHECKPOINT_FILE):
            cmonkey_run.run_from_checkpoint(CHECKPOINT_FILE)
        else: cmonkey_run.run()
