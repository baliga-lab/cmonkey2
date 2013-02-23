# vi: sw=4 ts=4 et:
"""human.py - cMonkey human specific module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import numpy as np
import scipy
import util
import thesaurus
import organism
import datamatrix as dm
import scoring
import microarray
import stringdb
import network as nw
import motif
import meme
import membership as memb
import set_enrichment as se
import cmonkey_run


CACHE_DIR = 'humancache'
CONTROLS_FILE = 'human_data/gbmTCGA/controls.csv'
RUG_FILE = 'human_data/gbmTCGA/rug.csv'

THESAURUS_FILE = 'human_data/gbmTCGA/synonymThesaurus.csv.gz'
PROM_SEQFILE = 'human_data/gbmTCGA/promoterSeqs_Homo_sapiens_gs.csv.gz'
P3UTR_SEQFILE = 'human_data/gbmTCGA/p3utrSeqs_Homo_sapiens_gs.csv.gz'


RUG_PROPS = ['GBM','CONTROL']
NUM_CLUSTERS = 134
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 3000

#SEQUENCE_TYPES = ['promoter', 'p3utr']
#SEARCH_DISTANCES = {'promoter': (0, 700), 'p3utr': (0, 831)}
#SCAN_DISTANCES = {'promoter': (0, 700), 'p3utr': (0, 831)}
#SEQ_FILENAMES = {'promoter': PROM_SEQFILE, 'p3utr': P3UTR_SEQFILE}
SEQUENCE_TYPES = ['upstream', 'p3utr']
SEARCH_DISTANCES = {'upstream': (0, 700), 'p3utr': (0, 831)}
SCAN_DISTANCES = {'upstream': (0, 700), 'p3utr': (0, 831)}
SEQ_FILENAMES = {'upstream': PROM_SEQFILE, 'p3utr': P3UTR_SEQFILE}

MAX_MOTIF_WIDTH = 12
SELECT_ROWS = True

# configure the function setup here
ADD_SET_ENRICHMENT = False
ADD_MEME = False
ADD_WEEDER = True

# scoring-specific
USE_SET_TYPES = ['pita'] # ['target_scan']
WEEDER_SEQ_TYPE = 'upstream'
MOTIF_START_ITERATION = 600
MOTIF_UPDATE_INTERVAL = 10
MOTIF_COMPUTE_INTERVAL = 100


def genes_per_group_proportional(num_genes_total, num_per_group):
    """takes the total number of genes and a dictionary containing
    the number of members for each group and returns a dictionary that
    distributes the number of genes proportionally to each group"""
    result = {}
    num_group_elems = sum(num_per_group.values())
    groups = num_per_group.keys()
    for index in range(len(groups)):
        group = groups[index]
        if index == len(groups) - 1:
            result[group] = num_genes_total - sum(result.values())
        else:
            result[group] = int(float(num_genes_total) *
                                float(num_per_group[group]) /
                                float(num_group_elems))
    return result


def genes_per_group_nonproportional(num_genes_total, groups):
    """distributes the number of genes evenly to each group given"""
    result = {}
    partition = int(float(num_genes_total) / float(len(groups)))
    for index in range(len(groups)):
        if index == len(groups) - 1:
            result[groups[index]] = num_genes_total - sum(result.values())
        else:
            result[groups[index]] = partition
    return result


def select_probes(matrix, num_genes_total, column_groups, proportional=True):
    """select probes proportional, column_groups is a map from a group
    label to column indexes in the matrix"""
    def coeff_var(row_values):
        """computes the coefficient of variation"""
        sigma = util.r_stddev(row_values)
        mu = util.mean(row_values)
        return sigma / mu

    num_per_group = {group: len(indexes)
                     for group, indexes in column_groups.items()}
    if proportional:
        per_group = genes_per_group_proportional(num_genes_total,
                                                 num_per_group)
    else:
        per_group = genes_per_group_proportional(num_genes_total,
                                                 column_groups.keys())

    cvrows = []
    mvalues = matrix.values
    nrows = matrix.num_rows
    for group, col_indexes in column_groups.items():
        group_cvs = []
        for row in xrange(nrows):
            row_values = [mvalues[row][col] for col in col_indexes]
            group_cvs.append(coeff_var(row_values))
        cvrows += [group_cvs.index(value)
                   for value in
                   sorted(group_cvs, reverse=True)][:per_group[group]]
    return sorted(list(set(cvrows)))


def intensities_to_ratios(matrix, controls, column_groups):
    """turn intensities into ratios
    Warning: input matrix is modified !!!"""
    colnames = matrix.column_names
    control_indexes = [colnames.index(control)
                       for control in controls
                       if control in colnames]
    mvalues = matrix.values
    nrows = matrix.num_rows
    for group_columns in column_groups.values():
        group_controls = [index for index in control_indexes
                          if index in group_columns]
        means = []
        for row in xrange(nrows):
            values = [float(mvalues[row][col]) for col in group_controls]
            means.append(sum(values) / float(len(values)))

        for col in group_columns:
            for row in range(matrix.num_rows):
                mvalues[row][col] /= means[row]

        center_scale_filter(matrix, group_columns, group_controls)
    return matrix


def center_scale_filter(matrix, group_columns, group_controls):
    """center the values of each row around their median and scale
    by their standard deviation. This is a specialized version"""
    mvalues = matrix.values
    centers = [scipy.median([mvalues[row][col]
                             for col in group_controls])
               for row in range(matrix.num_rows)]
    scale_factors = [util.r_stddev([mvalues[row][col]
                                    for col in group_columns])
                     for row in range(matrix.num_rows)]
    for row in range(matrix.num_rows):
        for col in group_columns:
            mvalues[row][col] -= centers[row]
            mvalues[row][col] /= scale_factors[row]
    return matrix


def read_controls():
    """reads the controls file"""
    with open(CONTROLS_FILE) as infile:
        return [line.strip() for line in infile.readlines()]


def read_rug(pred):
    """reads the rug file"""
    infile = util.read_dfile(RUG_FILE, sep=',', has_header=False)
    return list(set([row[0] for row in infile.lines if pred(row)]))


def read_matrix(filename):
    """reads the data matrix from a file"""
    controls = read_controls()
    rug = read_rug(lambda row: row[1] in RUG_PROPS)
    columns_to_use = list(set(rug + controls))

    # pass the column filter as the first filter to the DataMatrixFactory,
    # so normalization will be applied to the submatrix
    matrix_factory = dm.DataMatrixFactory([
            lambda matrix: matrix.submatrix_by_name(
                column_names=columns_to_use)])
    infile = util.read_dfile(filename, sep='\t', has_header=True, quote="\"")
    matrix = matrix_factory.create_from(infile)

    column_groups = {1: range(matrix.num_columns)}
    if SELECT_ROWS:
        #select_rows = select_probes(matrix, 2000, column_groups)
        #matrix = matrix.submatrix_by_rows(select_rows)
        matrix = matrix.submatrix_by_rows(range(2010))
    return intensities_to_ratios(matrix, controls, column_groups)


class RembrandtCMonkeyRun(cmonkey_run.CMonkeyRun):

    def __init__(self, organism_code, ratio_matrix, num_clusters):
        cmonkey_run.CMonkeyRun.__init__(self, organism_code, ratio_matrix,
                                        num_clusters=num_clusters)
        self.__organism = None
        self['cache_dir'] = CACHE_DIR
        self['sequence_types'] = SEQUENCE_TYPES
        self['search_distances'] = SEARCH_DISTANCES
        self['scan_distances'] = SCAN_DISTANCES
        self['num_iterations'] = NUM_ITERATIONS

    def organism(self):
        if self.__organism == None:
            self.__organism = self.make_hsa()
        return self.__organism

    def make_hsa(self):
        """returns a human organism object"""
        nw_factories = [stringdb.get_network_factory3('human_data/gbmTCGA/string_geneSymbol.csv', weight=1.0)]
        return organism.GenericOrganism('hsa', THESAURUS_FILE, nw_factories,
                                        seq_filenames=SEQ_FILENAMES,
                                        search_distances=SEARCH_DISTANCES,
                                        scan_distances=SCAN_DISTANCES)

    def meme_suite(self, seqtype):
        """upstream meme suite"""
        background_file = meme.global_background_file(
            self.organism(), self.ratio_matrix.row_names, seqtype, use_revcomp=True)
        return meme.MemeSuite430(max_width=MAX_MOTIF_WIDTH,
                                 background_file=background_file)

    def make_meme_scoring(self, seqtype, meme_suite, sequence_filters):
        motif_scaling_fun = scoring.get_default_motif_scaling(self['num_iterations'])
        return motif.MemeScoringFunction(
            self.organism(), self.membership(),
            self.ratio_matrix, meme_suite,
            seqtype=seqtype,
            sequence_filters=sequence_filters,
            scaling_func=motif_scaling_fun,
            num_motif_func=motif.default_nmotif_fun,
            update_in_iteration=scoring.schedule(MOTIF_START_ITERATION, MOTIF_UPDATE_INTERVAL),
            motif_in_iteration=scoring.schedule(MOTIF_START_ITERATION, MOTIF_COMPUTE_INTERVAL),
            config_params=self.config_params)

    def make_weeder_scoring(self, seqtype, meme_suite, sequence_filters):
        motif_scaling_fun = scoring.get_default_motif_scaling(self['num_iterations'])
        return motif.WeederScoringFunction(
            self.organism(), self.membership(), self.ratio_matrix,
            meme_suite,
            seqtype=seqtype,
            sequence_filters=sequence_filters,
            scaling_func=motif_scaling_fun,
            num_motif_func=motif.default_nmotif_fun,
            update_in_iteration=scoring.schedule(MOTIF_START_ITERATION, MOTIF_UPDATE_INTERVAL),
            motif_in_iteration=scoring.schedule(MOTIF_START_ITERATION, MOTIF_COMPUTE_INTERVAL),
            config_params=self.config_params)

    def make_network_scoring(self, scaling_fun):
        return nw.ScoringFunction(self.organism(),
                                  self.membership(),
                                  self.ratio_matrix,
                                  scaling_func=scaling_fun,
                                  run_in_iteration=scoring.schedule(1, 7),
                                  config_params=self.config_params)

    def make_set_scoring(self, scaling_fun, use_set_types):
        set_types = []
        if 'pita' in use_set_types:
            set_types.append(se.SetType.read_csv('pita', 'human_data/gbmTCGA/pita_miRNA_sets_geneSymbol.csv'))
        if 'target_scan' in use_set_types:
            set_types.append(se.SetType.read_csv('target_scan',
                                                 'human_data/gbmTCGA/targetScan_miRNA_sets_geneSymbol.csv'))
        return se.ScoringFunction(self.membership(), self.ratio_matrix,
                                  set_types, scaling_fun,
                                  scoring.schedule(1, 7),
                                  config_params=self.config_params)

    def make_row_scoring(self):
        """returns the row scoring function"""
        # we need to remove the location from the sequence when selecting for
        # individual clusters
        sequence_filters = [lambda seqs, feature_ids: {key: seqs[key][1] for key in seqs.keys()}]
        network_scaling_fun = scoring.get_default_network_scaling(self['num_iterations'])
        meme_scoring = None
        weeder_scoring = None

        row_scoring = microarray.RowScoringFunction(
            self.membership(), self.ratio_matrix,
            lambda iteration: ROW_WEIGHT,
            config_params=self.config_params)
        scoring_funcs = [row_scoring,
                         self.make_network_scoring(network_scaling_fun)]

        if ADD_SET_ENRICHMENT:
            scoring_funcs.append(self.make_set_scoring(network_scaling_fun,
                                                       USE_SET_TYPES))

        if ADD_MEME:
            meme_suite_meme = self.meme_suite('upstream')
            meme_scoring = self.make_meme_scoring('upstream',
                                                  meme_suite_meme,
                                                  sequence_filters + [motif.get_remove_low_complexity_filter(meme_suite_meme)])

        if ADD_WEEDER:
            meme_suite_weeder = self.meme_suite(WEEDER_SEQ_TYPE)
            weeder_scoring = self.make_weeder_scoring(WEEDER_SEQ_TYPE,
                                                      meme_suite_weeder,
                                                      sequence_filters + [motif.get_remove_low_complexity_filter(meme_suite_weeder)])
        if ADD_MEME and ADD_WEEDER:
            scoring_funcs.append(scoring.ScoringFunctionCombiner(
                    self.membership(),
                    [meme_scoring, weeder_scoring],
                    scaling_func=lambda iteration: 0.5,
                    config_params=self.config_params))
        else:
            if ADD_MEME:
                scoring_funcs.append(meme_scoring)
            if ADD_WEEDER:
                scoring_funcs.append(weeder_scoring)

        return scoring.ScoringFunctionCombiner(self.membership(), scoring_funcs,
                                               config_params=self.config_params)


if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011-2012, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) < 2:
        print('Usage: ./gbmTCGA.sh <ratio-file> [checkpoint-file]')
    else:
        if len(sys.argv) > 2:
            CHECKPOINT_FILE = sys.argv[2]
        cmonkey_run = RembrandtCMonkeyRun('hsa', read_matrix(sys.argv[1]), NUM_CLUSTERS)
        # Turns of multiprocessing for debugging purposes
        #cmonkey_run['multiprocessing'] = False
        cmonkey_run.run()
