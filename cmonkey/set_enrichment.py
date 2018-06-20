# vi: sw=4 ts=4 et:
"""set_enrichment.py - cMonkey set_enrichment scoring.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import math
import os
import json
import logging
import numpy as np
import multiprocessing as mp
from collections import defaultdict

import cmonkey.util as util
import cmonkey.scoring as scoring
import cmonkey.datamatrix as dm

# Python2/Python3 compatibility
try:
    xrange
except NameError:
    xrange = range


class DiscreteEnrichmentSet:
    def __init__(self, genes):
        """instance creation"""
        self.__genes = genes
        self.cutoff = "discrete"

    def genes(self):
        return self.__genes

    def genes_above_cutoff(self):
        """returns the genes that have a weight above the cutoff"""
        return self.genes()

    def __repr__(self):
        return ("Enrichment set: Cutoff = discrete, # genes: %d # above cutoff: %d" %
                (len(self.__genes), len(self.__genes)))


class CutoffEnrichmentSet:
    """Enrichment set representation constructed with a cutoff"""
    def __init__(self, cutoff, elems):
        """instance creation"""
        self.elems = elems  # pairs of gene, weight
        self.cutoff = cutoff
        self.__genes = None

    def genes(self):
        """returns all genes"""
        if self.__genes is None:
            self.__genes = {elem[0] for elem in self.elems}
        return self.__genes

    def genes_above_cutoff(self):
        """returns the genes that have a weight above the cutoff"""
        return {elem[0] for elem in self.elems if elem[1] >= self.cutoff}

    def __repr__(self):
        return ("Enrichment set: Cutoff = %s, # genes: %d # above cutoff: %d" %
                (self.cutoff, len(self.elems), len(self.genes_above_cutoff())))


class SetType:
    """Set type representation. This is just a grouping from name to a number of sets
    and providing access to all contained genes"""

    def __init__(self, name, sets, weight):
        """instance creation"""
        self.name = name
        self.sets = sets
        self.weight = weight
        self.__genes = None

    def genes(self):
        """All genes contained in the sets"""
        if self.__genes is None:
            all_genes = [enrichment_set.genes() for enrichment_set in self.sets.values()]
            self.__genes = set.union(*all_genes)
        return self.__genes

    def __repr__(self):
        """string representation"""
        result = "SetType['%s'] = {\n" % self.name
        for key, value in self.sets.items():
            result += "%s -> %s\n" % (key, value)
        result += "}"
        return result


# global variables that are shared by the child processe
# when running in multiprocessing. This is an attempt to avoid
# that the scoring function slows itself down by using to much
# memory
SET_MATRIX = None
SET_MEMBERSHIP = None
SET_SET_TYPE = None
SET_SYNONYMS = None
CANONICAL_ROWNAMES = None
CANONICAL_ROW_INDEXES = None


def read_set_types(config_params, thesaurus, input_genes):
    """Reads sets from a JSON file. We also ensure that genes
    are stored in canonical form in the set, so that set operations based on
    gene names will succeed"""
    result = []
    set_types = config_params['SetEnrichment']['set_types'].split(',')

    for set_type in set_types:
        setfile = config_params['SetEnrichment-%s' % set_type]['set_file']
        weight = float(config_params['SetEnrichment-%s' % set_type]['weight'])
        if setfile.endswith('csv'):
            with open(setfile) as infile:
                sets = read_sets_csv(infile)
        else:
            with open(setfile) as infile:
                sets = json.load(infile)
        sets = process_sets(sets, thesaurus, input_genes)
        result.append(SetType(set_type, sets, weight))
    return result


def process_sets(input_sets, thesaurus, input_genes):
    """Reusable function that maps a dictionary {name: [genes]}
    to a set of DiscretenEnrichmentSet objects"""
    sets = {}
    genes_thrown_out = 0
    sets_thrown_out = 0
    input_genes = {thesaurus[gene] for gene in input_genes if gene in thesaurus}

    for setname, genes in input_sets.items():
        filtered = map(str, filter(lambda g: g in thesaurus, genes))
        canonic_genes = {thesaurus[gene] for gene in filtered}
        genes_thrown_out += len(genes) - len(canonic_genes)

        # check whether the genes are found in the ratios
        # and ignore
        canonic_genes = {gene for gene in canonic_genes if gene in input_genes}
        if len(canonic_genes) > 0:
            sets[setname] = DiscreteEnrichmentSet(canonic_genes)
        else:
            sets_thrown_out += 1
    logging.info("SET_ENRICHMENT REMOVED %d GENES FROM INPUT", genes_thrown_out)
    logging.info("SET_ENRICHMENT REMOVED %d SETS FROM INPUT", sets_thrown_out)
    return sets


def read_sets_csv(infile, sep1=',', sep2=';'):
    """Reads sets from a CSV file
    We support 2-column and 3 column formats:

    2 column files have the format

    <set name><separator><gene>

    3 column files have the format

    <set name><separator><gene><weight>
    """
    line1 = infile.readline().strip().split(sep1)
    if len(line1) == 2:
        sets = defaultdict(list)
        for gene in line1[1].split(sep2):
            sets[line1[0]].append(gene)

        for line in infile:
            row = line.strip().split(sep1)
            for gene in row[1].split(sep2):
                sets[row[0]].append(gene)
        return sets
    else:
        raise Exception("3 column set files not supported yet")


class ScoringFunction(scoring.ScoringFunctionBase):
    """Set enrichment scoring function"""
    def __init__(self, function_id, cmrun):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, function_id, cmrun)
        self.__set_types = read_set_types(self.config_params, self.organism.thesaurus(),
                                          self.ratios.row_names)
        self.run_log = scoring.RunLog(function_id, cmrun.dbsession(), self.config_params)

    def bonferroni_cutoff(self):
        """Bonferroni cutoff value"""
        return 0.05 / float(self.num_clusters())

    def do_compute(self, iteration_result, ref_matrix):
        """compute method
        Note: will return None if not computed yet and the result of a previous
        scoring if the function is not supposed to actually run in this iteration
        """
        global SET_MATRIX, SET_MEMBERSHIP, SET_SET_TYPE, SET_SYNONYMS, CANONICAL_ROWNAMES, CANONICAL_ROW_INDEXES
        logging.info("Compute scores for set enrichment...")
        start_time = util.current_millis()
        matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                               self.gene_names())
        use_multiprocessing = self.config_params[scoring.KEY_MULTIPROCESSING]
        SET_MATRIX = self.ratios
        SET_MEMBERSHIP = self.membership
        SET_SYNONYMS = self.organism.thesaurus()

        if CANONICAL_ROWNAMES is None:
            CANONICAL_ROWNAMES = set(map(lambda n: SET_SYNONYMS[n] if n in SET_SYNONYMS else n,
                                         self.ratios.row_names))

        if CANONICAL_ROW_INDEXES is None:
            CANONICAL_ROW_INDEXES = {}
            for index, row in enumerate(self.ratios.row_names):
                if row in SET_SYNONYMS:
                    CANONICAL_ROW_INDEXES[SET_SYNONYMS[row]] = index
                else:
                    CANONICAL_ROW_INDEXES[row] = index

        ref_min_score = np.nanpercentile(ref_matrix.values, 10.0)
        logging.info('REF_MIN_SCORE: %f', ref_min_score)

        set_filepath = os.path.join(self.config_params['output_dir'],
                                    'setEnrichment_set.csv')
        pval_filepath = os.path.join(self.config_params['output_dir'],
                                     'setEnrichment_pvalue.csv')

        for set_type in self.__set_types:
            SET_SET_TYPE = set_type
            logging.info("PROCESSING SET TYPE '%s'", set_type.name)
            start1 = util.current_millis()
            cutoff = self.bonferroni_cutoff()
            if use_multiprocessing:
                with util.get_mp_pool(self.config_params) as pool:
                    results = pool.map(compute_cluster_score,
                                       [(cluster, cutoff, ref_min_score)
                                        for cluster in xrange(1, self.num_clusters() + 1)])
            else:
                results = []
                for cluster in xrange(1, self.num_clusters() + 1):
                    results.append(compute_cluster_score((cluster, cutoff, ref_min_score)))

            elapsed1 = util.current_millis() - start1
            logging.info("ENRICHMENT SCORES COMPUTED in %f s, STORING...",
                         elapsed1 / 1000.0)

            if not os.path.exists(set_filepath):
                setFile = open(set_filepath, 'w')
                setFile.write(',' + ','.join([str(i) for i in xrange(1, self.num_clusters() + 1)]))
                pvFile = open(pval_filepath, 'w')
                pvFile.write(',' + ','.join([str(i) for i in xrange(1, self.num_clusters() + 1)]))
            else:
                setFile = open(set_filepath, 'a')
                pvFile = open(pval_filepath, 'a')

            minSets = []
            pValues = []
            for cluster in xrange(1, self.num_clusters() + 1):
                # store the best enriched set determined
                scores, min_set, min_pvalue = results[cluster - 1]
                minSets.append(min_set)
                pValues.append(min_pvalue)

                for row in xrange(len(self.gene_names())):
                    matrix.values[row][cluster - 1] += scores[row] * set_type.weight
            setFile.write('\n'+str(iteration_result['iteration'])+','+','.join([str(i) for i in minSets]))
            pvFile.write('\n'+str(iteration_result['iteration'])+','+','.join([str(i) for i in pValues]))
            setFile.close()
            pvFile.close()

        logging.info("SET ENRICHMENT FINISHED IN %f s.\n",
                     (util.current_millis() - start_time) / 1000.0)
        # cleanup
        SET_SET_TYPE = None
        SET_MATRIX = None
        SET_MEMBERSHIP = None
        SET_SYNONYMS = None

        return matrix

    def run_logs(self):
        """return the run logs"""
        return [self.run_log]


def compute_cluster_score(args):
    """Computes the cluster score for a given set type"""
    global SET_MATRIX, SET_MEMBERSHIP, SET_SET_TYPE, SET_SYNONYMS, CANONICAL_ROWNAMES, CANONICAL_ROW_INDEXES
    cluster, cutoff, ref_min_score = args
    return compute_cluster_score_plain(cluster, cutoff, ref_min_score, SET_MATRIX, SET_MEMBERSHIP,
                                       SET_SET_TYPE, SET_SYNONYMS, CANONICAL_ROWNAMES, CANONICAL_ROW_INDEXES)

def compute_cluster_score_plain(cluster, cutoff, ref_min_score, SET_MATRIX, SET_MEMBERSHIP, SET_SET_TYPE,
                                SET_SYNONYMS, CANONICAL_ROWNAMES, CANONICAL_ROW_INDEXES):
    """This version is the real implementation, that can be tested without using global variables"""
    set_type = SET_SET_TYPE
    matrix = SET_MATRIX
    cluster_rows = set()
    for gene in SET_MEMBERSHIP.rows_for_cluster(cluster):
        if gene in SET_SYNONYMS:
            cluster_rows.add(SET_SYNONYMS[gene])
        else:
            cluster_rows.add(gene)
    set_type_genes = set_type.genes()

    cluster_genes = {gene for gene in cluster_rows if gene in set_type_genes}
    overlap_sizes = []
    set_sizes = []
    set_names = []
    for set_name in sorted(set_type.sets.keys()):
        eset = set_type.sets[set_name]
        set_genes = eset.genes_above_cutoff()
        intersect = len(cluster_genes.intersection(set_genes))
        if intersect > 0:
            set_names.append(set_name)
            set_sizes.append(len(set_genes))
            overlap_sizes.append(intersect)

    num_sets = len(overlap_sizes)
    scores = np.zeros(matrix.num_rows)
    if not num_sets==0:
        num_genes = len(set_type_genes)
        phyper_n = list(np.array([num_genes] * num_sets) - np.array(set_sizes))
        phyper_k = [len(cluster_genes)] * num_sets

        enrichment_pvalues = np.array(util.phyper(overlap_sizes, set_sizes, phyper_n, phyper_k))
        min_pvalue = enrichment_pvalues[np.isfinite(enrichment_pvalues)].min()
        min_index = np.where(enrichment_pvalues == min_pvalue)[0][0]
        min_set = set_names[min_index]
        min_set_overlap = overlap_sizes[min_index]

        min_genes = set_type.sets[min_set].genes()
        # ensure all row names are in canonical form
        min_genes = [gene for gene in min_genes if gene in CANONICAL_ROWNAMES]
        min_indexes = [CANONICAL_ROW_INDEXES[gene] for gene in min_genes]

        if set_type.sets[min_set].cutoff == 'discrete':
            overlap_genes = cluster_genes.intersection(set(min_genes))
            overlap_indexes = [CANONICAL_ROW_INDEXES[gene] for gene in overlap_genes]
            scores[min_indexes] = 0.5
            scores[overlap_indexes] = 1.0
        else:
            # scaling
            for index in min_indexes:
                scores[index] = set_type.sets[min_set].weights[index]
            min_weight = min_set_weights.min()
            max_weight = min_set_weights.max()
            weight_range = max_weight - min_weight
            scores[min_indexes] -= min_weight
            scores[min_indexes] /= weight_range

        if min_pvalue <= cutoff:
            dampened_pvalue = 1
        else:
            dampened_pvalue = math.log10(min_pvalue) / math.log10(cutoff)

        scores[scores != 0.0] = scores[scores != 0.0] / dampened_pvalue
        scores *= ref_min_score
    else:
        min_set = 'NA'
        min_pvalue = np.nan

    return scores, min_set, min_pvalue
