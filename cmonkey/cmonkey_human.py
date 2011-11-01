"""cmonkey_human.py - cMonkey human module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import os
import logging
import datamatrix as dm
import util
import human
import scipy.cluster.vq as clvq
import membership as memb
import microarray
import stringdb
import network as nw
import motif
import meme

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
CMONKEY_VERSION = '4.0'
CACHE_DIR = 'humancache'
CONTROLS_FILE = 'human_data/controls.csv'
RUG_FILE = 'human_data/rug.csv'

PROM_SEQFILE = 'human_data/promoterSeqs_set3pUTR_Final.csv.gz'
P3UTR_SEQFILE = 'human_data/p3utrSeqs_set3pUTR_Final.csv.gz'
THESAURUS_FILE = 'human_data/synonymThesaurus.csv.gz'

RUG_PROPS = ['MIXED', 'ASTROCYTOMA', 'GBM', 'OLIGODENDROGLIOMA']
NUM_CLUSTERS = 133
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 2000

def run_cmonkey():
    """init of the cMonkey system"""
    logging.basicConfig(format=LOG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    if not os.path.exists(CACHE_DIR):
        os.mkdir(CACHE_DIR)
    matrix = read_matrix(sys.argv[1])
    logging.info("Filtered input matrix has %d rows and %d columns:",
                 matrix.num_rows(), matrix.num_columns())
    membership = make_membership(matrix)
    logging.info("Memberships seeded")
    organism = make_organism()
    logging.info("Organism object created")
    gene_scoring_funcs = make_gene_scoring_funcs(organism, membership, matrix)
    cond_scoring = microarray.ColumnScoringFunction(membership, matrix)
    logging.info("Scoring functions created")

    for iteration in range(NUM_ITERATIONS):
        iterate(membership, matrix, gene_scoring_funcs, cond_scoring,
                iteration, NUM_ITERATIONS)
    logging.info("Done !!!!")


def iterate(membership, matrix, gene_scoring_funcs, cond_scoring_func,
            iteration, num_iterations):
    """one iteration of the algorithm"""
    logging.info("Iteration # %d", iteration)
    result_matrices = []
    score_weights = []
    logging.info("calculating row scores...")
    for score_func in gene_scoring_funcs:
        row_scores = score_func.compute(iteration)
        if row_scores != None:
            result_matrices.append(row_scores)
            score_weights.append(score_func.weight(iteration))
    logging.info("row scores done.")
    logging.info("calculating column scores...")
    cscores = cond_scoring_func.compute(iteration)
    logging.info("column scores done.")

    if len(result_matrices) > 0:  # should be 0
        logging.info("quantile normalize...")
        result_matrices = dm.quantile_normalize_scores(result_matrices,
                                                       score_weights)
        logging.info("quantile normalize done.")

    rscores = result_matrices[0] * gene_scoring_funcs[0].weight(iteration)
    for index in range(1, len(result_matrices)):
        rscores = rscores + (result_matrices[index] *
                             gene_scoring_funcs[index].weight(iteration))

    membership.update(matrix, rscores, cscores, iteration, num_iterations)


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
    select_rows = human.select_probes(matrix, 2000, column_groups)
    matrix = matrix.submatrix_by_rows(select_rows)
    return human.intensities_to_ratios(matrix, controls, column_groups)


def make_membership(matrix):
    """returns a seeded membership"""
    return memb.ClusterMembership.create(
        matrix.sorted_by_row_name(),
        memb.make_kmeans_row_seeder(NUM_CLUSTERS),
        microarray.seed_column_members,
        num_clusters=NUM_CLUSTERS)

def make_organism():
    """returns a human organism object"""
    nw_factories = [stringdb.get_network_factory3('human_data/string.csv')]
    organism = human.Human(PROM_SEQFILE, P3UTR_SEQFILE, THESAURUS_FILE,
                           nw_factories)
    return organism


def make_gene_scoring_funcs(organism, membership, matrix):
    """setup the gene-related scoring functions here
    each object in this array supports the method
    compute(organism, membership, matrix) and returns
    a DataMatrix(genes x cluster)
    """
    row_scoring = microarray.RowScoringFunction(membership, matrix,
                                                lambda iteration: ROW_WEIGHT)

    sequence_filters = []
    meme_suite = meme.MemeSuite430()
    motif_scoring = motif.ScoringFunction(organism,
                                          membership,
                                          matrix,
                                          meme_suite,
                                          sequence_filters=sequence_filters,
                                          pvalue_filter=motif.make_min_value_filter(-20.0),
                                          weight_func=lambda iteration: 0.0,
                                          interval=1)

    network_scoring = nw.ScoringFunction(organism, membership, matrix,
                                         lambda iteration: 0.0, 7)
    return [row_scoring, motif_scoring, network_scoring]

if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    print('(Human version)')
    if len(sys.argv) <= 1:
        print('Usage: ./run_cmonkey_human.sh <ratio-file>')
    else:
        run_cmonkey()
