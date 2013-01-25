# vi: sw=4 ts=4 et:
"""leishmania.py - Mouse/Leishmania configuration

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
import cmonkey
import meme
import motif
import stringdb
import organism
import network as nw
import cmonkey_run

CMONKEY_VERSION = '4.0'
CHECKPOINT_INTERVAL = 3
CHECKPOINT_FILE = None

NUM_ITERATIONS = 2000
CACHE_DIR = 'leishmania_cache'
NUM_CLUSTERS = 267
ROW_WEIGHT = 6.0
MAX_MOTIF_WIDTH = 12

THESAURUS_FILE = 'leishmania_data/synonymThesaurus.csv.gz'

SEQUENCE_TYPES = ['upstream', 'p3utr']
SEARCH_DISTANCES = {'upstream': (0, 700), 'p3utr': (0, 831)}
SCAN_DISTANCES = {'upstream': (0, 700), 'p3utr': (0, 831)}
PROM_SEQFILE = 'leishmania_data/promoterSeqs_Mus_musculus_RA.csv'
P3UTR_SEQFILE = 'leishmania_data/p3utrSeqs_Mus_musculus_RA.csv.gz'
SEQ_FILENAMES = {'upstream': PROM_SEQFILE, 'p3utr': P3UTR_SEQFILE}


class LeishmaniaCMonkeyRun(cmonkey_run.CMonkeyRun):
    def __init__(self, organism_code, ratio_matrix, num_clusters):
        cmonkey_run.CMonkeyRun.__init__(self, organism_code, ratio_matrix, num_clusters)
        self.__organism = None
        self['cache_dir'] = CACHE_DIR
        self['sequence_types'] = SEQUENCE_TYPES
        self['search_distances'] = SEARCH_DISTANCES
        self['scan_distances'] = SCAN_DISTANCES

    def organism(self):
        if self.__organism == None:
            self.__organism = self.make_mmu()
        return self.__organism

    def make_mmu(self):
        """returns a configured organism object"""
        nw_factories = [stringdb.get_network_factory2(
                'mmu',
                'leishmania_data/mouse_links_preprocessed.csv', weight=1.0, sep=';')]
        return organism.GenericOrganism('mmu', THESAURUS_FILE, nw_factories,
                                        seq_filenames=SEQ_FILENAMES,
                                        search_distances=SEARCH_DISTANCES,
                                        scan_distances=SCAN_DISTANCES)

    def make_row_scoring(self):
        """returns the row scoring function"""
        row_scoring = microarray.RowScoringFunction(
            self.membership(), self.ratio_matrix,
            lambda iteration: ROW_WEIGHT,
            config_params=self.config_params)

        # motifing
        sequence_filters = []
        background_file_prom = meme.global_background_file(
            self.organism(), self.ratio_matrix.row_names(), 'upstream',
            use_revcomp=True)
        background_file_p3utr = meme.global_background_file(
            self.organism(), self.ratio_matrix.row_names(), 'p3utr',
            use_revcomp=True)
        meme_suite_prom = meme.MemeSuite430(
            max_width=MAX_MOTIF_WIDTH,
            background_file=background_file_prom)
        meme_suite_p3utr = meme.MemeSuite430(
            max_width=MAX_MOTIF_WIDTH,
            background_file=background_file_p3utr)

        motif_scoring = motif.MemeScoringFunction(
            self.organism(),
            self.membership(),
            self.ratio_matrix,
            meme_suite_prom,
            seqtype='upstream',
            sequence_filters=sequence_filters,
            scaling_func=lambda iteration: 0.0,
            run_in_iteration=scoring.default_motif_iterations,
            config_params=self.config_params)

        weeder_scoring = motif.WeederScoringFunction(
            self.organism(), self.membership(), self.ratio_matrix,
            meme_suite_p3utr, 'p3utr',
            scaling_func=lambda iteration: 0.0,
            run_in_iteration=scoring.default_motif_iterations,
            config_params=self.config_params)

        motif_combiner = scoring.ScoringFunctionCombiner(
            self.membership(), [motif_scoring, weeder_scoring],
            scaling_func=lambda iteration: 0.5,
            config_params=self.config_params)

        network_scoring = nw.ScoringFunction(self.organism(),
                                             self.membership(),
                                             self.ratio_matrix,
                                             lambda iteration: 0.0,
                                             scoring.default_network_iterations,
                                             config_params=self.config_params)

        return scoring.ScoringFunctionCombiner(
            self.membership(), [row_scoring, network_scoring, motif_combiner],
            config_params=self.config_params)


if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) < 2:
        print('Usage: ./leishmania.sh <ratio-file> [checkpoint-file]')
    else:
        if len(sys.argv) > 2:
            CHECKPOINT_FILE = sys.argv[2]

        matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
        infile = util.read_dfile(sys.argv[1], has_header=True, quote='\"')
        matrix = matrix_factory.create_from(infile)
        cmonkey_run = LeishmaniaCMonkeyRun('mmu', matrix, NUM_CLUSTERS)
        cmonkey_run.run()
