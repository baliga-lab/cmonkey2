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
CONTROLS_FILE = 'human_data/controls.csv'
RUG_FILE = 'human_data/rug.csv'

THESAURUS_FILE = 'human_data/synonymThesaurus.csv.gz'

RUG_PROPS = ['MIXED', 'ASTROCYTOMA', 'GBM', 'OLIGODENDROGLIOMA']
NUM_CLUSTERS = 720  # 133
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 2000

SEQUENCE_TYPES = ['promoter', 'p3utr']
SEARCH_DISTANCES = {'promoter': (0, 700), 'p3utr': (0, 831)}
SCAN_DISTANCES = {'promoter': (0, 700), 'p3utr': (0, 831)}
PROM_SEQFILE = 'human_data/promoterSeqs_set3pUTR_Final.csv.gz'
P3UTR_SEQFILE = 'human_data/p3utrSeqs_set3pUTR_Final.csv.gz'
SEQ_FILENAMES = {'promoter': PROM_SEQFILE, 'p3utr': P3UTR_SEQFILE}
MAX_MOTIF_WIDTH = 12


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
    for group, col_indexes in column_groups.items():
        group_cvs = []
        for row in range(matrix.num_rows()):
            row_values = [matrix[row][col] for col in col_indexes]
            group_cvs.append(coeff_var(row_values))
        cvrows += [group_cvs.index(value)
                   for value in
                   sorted(group_cvs, reverse=True)][:per_group[group]]
    return sorted(list(set(cvrows)))


def intensities_to_ratios(matrix, controls, column_groups):
    """turn intensities into ratios
    Warning: input matrix is modified !!!"""
    colnames = matrix.column_names()
    control_indexes = [np.where(colnames == control)[0][0]
                       for control in controls
                       if control in colnames]
    for group_columns in column_groups.values():
        group_controls = [index for index in control_indexes
                          if index in group_columns]
        means = []
        for row in range(matrix.num_rows()):
            values = [float(matrix[row][col]) for col in group_controls]
            means.append(sum(values) / float(len(values)))

        for col in group_columns:
            for row in range(matrix.num_rows()):
                matrix[row][col] /= means[row]

        center_scale_filter(matrix, group_columns, group_controls)
    return matrix


def center_scale_filter(matrix, group_columns, group_controls):
    """center the values of each row around their median and scale
    by their standard deviation. This is a specialized version"""
    centers = [scipy.median([matrix[row][col]
                             for col in group_controls])
               for row in range(matrix.num_rows())]
    scale_factors = [util.r_stddev([matrix[row][col]
                                    for col in group_columns])
                     for row in range(matrix.num_rows())]
    for row in range(matrix.num_rows()):
        for col in group_columns:
            matrix[row][col] -= centers[row]
            matrix[row][col] /= scale_factors[row]
    return matrix



def read_controls():
    """reads the controls file"""
    with open(CONTROLS_FILE) as infile:
        return [line.strip() for line in infile.readlines()]


def read_rug(pred):
    """reads the rug file"""
    infile = util.DelimitedFile.read(RUG_FILE, sep=',', has_header=False)
    return list(set([row[0] for row in infile.lines() if pred(row)]))


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
    infile = util.DelimitedFile.read(filename, sep=',', has_header=True,
                                     quote="\"")
    matrix = matrix_factory.create_from(infile)

    column_groups = {1: range(matrix.num_columns())}
    #select_rows = select_probes(matrix, 2000, column_groups)
    #matrix = matrix.submatrix_by_rows(select_rows)
    return intensities_to_ratios(matrix, controls, column_groups)


class RembrandtCMonkeyRun(cmonkey_run.CMonkeyRun):

    def __init__(self, organism_code, ratio_matrix, num_clusters):
        cmonkey_run.CMonkeyRun.__init__(self, organism_code, ratio_matrix, num_clusters)
        self.__organism = None
        self['cache_dir'] = CACHE_DIR
        self['sequence_types'] = SEQUENCE_TYPES
        self['search_distances'] = SEARCH_DISTANCES
        self['scan_distances'] = SCAN_DISTANCES

    def organism(self):
        if self.__organism == None:
            self.__organism = self.make_hsa()
        return self.__organism

    def make_hsa(self):
        """returns a human organism object"""
        nw_factories = [stringdb.get_network_factory3('human_data/string.csv', weight=1.0)]
        return organism.GenericOrganism('hsa', THESAURUS_FILE, nw_factories,
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
            self.organism(), self.ratio_matrix.row_names(), 'promoter',
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
            seqtype='promoter',
            sequence_filters=sequence_filters,
            pvalue_filter=motif.MinPValueFilter(-20.0),
            weight_func=lambda iteration: 0.0,
            run_in_iteration=scoring.default_motif_iterations,
            config_params=self.config_params)

        network_scoring = nw.ScoringFunction(self.organism(),
                                             self.membership(),
                                             self.ratio_matrix,
                                             lambda iteration: 0.0,
                                             scoring.default_network_iterations,
                                             config_params=self.config_params)

        weeder_scoring = motif.WeederScoringFunction(
            self.organism(), self.membership(), self.ratio_matrix,
            meme_suite_p3utr, 'p3utr',
            pvalue_filter=motif.MinPValueFilter(-20.0),
            weight_func=lambda iteration: 0.0,
            run_in_iteration=scoring.default_motif_iterations,
            config_params=self.config_params)

        pita = se.SetType.read_csv('pita', 'human_data/pita_miRNA_sets.csv')
        target_scan = se.SetType.read_csv(
            'target_scan', 'human_data/targetScan_miRNA_sets.csv')
        set_types = [pita, target_scan]
        set_enrichment_scoring = se.ScoringFunction(
            self.membership(),
            self.ratio_matrix,
            set_types,
            lambda iteration: 0.0, 7,
            config_params=self.config_params)

        motif_combiner = scoring.ScoringFunctionCombiner(
            self.membership(),
            [motif_scoring, weeder_scoring],
            weight_func=lambda iteration: 0.5)

        return scoring.ScoringFunctionCombiner(
            self.membership(),
            [row_scoring, motif_combiner, network_scoring,
             set_enrichment_scoring])


if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) < 2:
        print('Usage: ./rembrandt.sh <ratio-file> [checkpoint-file]')
    else:
        if len(sys.argv) > 2:
            CHECKPOINT_FILE = sys.argv[2]
        cmonkey_run = RembrandtCMonkeyRun('hsa', read_matrix(sys.argv[1]), NUM_CLUSTERS)
        cmonkey_run.run()
