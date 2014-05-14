# vi: sw=4 ts=4 et:
"""set_enrichment.py - cMonkey set_enrichment scoring.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import util
import math
import os
import json
import logging
import numpy as np

import scoring
import datamatrix as dm
import multiprocessing as mp


class EnrichmentSet:
    """Enrichment set representation"""
    def __init__(self, cutoff, elems):
        """instance creation"""
        self.elems = elems  # pairs of gene, weight
        self.cutoff = cutoff
        self.__genes_above_cutoff = None
        self.__genes = None

    def genes(self):
        if self.__genes is None:
            self.__genes = [elem[0] for elem in self.elems]
        return self.__genes

    def genes_above_cutoff(self):
        """returns the genes that have a weight above the cutoff"""
        if self.__genes_above_cutoff is None:
            if self.cutoff == 'discrete':
                self.__genes_above_cutoff = [elem[0] for elem in self.elems]
            else:
                self.__genes_above_cutoff = [elem[0] for elem in self.elems
                                             if elem[1] >= self.cutoff]            
        return self.__genes_above_cutoff

    def __repr__(self):
        return ("Enrichment set: Cutoff = %s, # genes: %d # above cutoff: %d" %
                (self.cutoff, len(self.elems), len(self.genes_above_cutoff())))


class SetType:
    """Set type representation. This is just a grouping from name to a number of sets
    and providing access to all contained genes"""

    def __init__(self, name, sets):
        """instance creation"""
        self.name = name
        self.sets = sets
        self.__genes = None

    def genes(self):
        """All genes contained in the sets"""
        if self.__genes is None:
            self.__genes = {elem[0] for enrichment_set in self.sets.values()
                            for elem in enrichment_set.elems}
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

def read_set_types(config_params):
    """Reads sets from a JSON file"""
    setfile = config_params['SetEnrichment']['set_file']
    with open(setfile) as infile:
        json_sets = json.load(infile)
        sets = {setname: EnrichmentSet('discrete', zip(genes, [1] * len(genes)))
                for setname, genes in json_sets.items()}
    return [SetType('default', sets)]


class ScoringFunction(scoring.ScoringFunctionBase):
    """Network scoring function"""
    def __init__(self, organism, membership, ratios, config_params=None):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, "SetEnrichment", organism, membership,
                                             ratios, config_params)
        self.__set_types = read_set_types(config_params)
        # stores (min_set, pvalue) pairs for each cluster and set type
        # for the last run of the function
        self.run_log = scoring.RunLog('set_enrichment', config_params)

    def bonferroni_cutoff(self):
        """Bonferroni cutoff value"""
        return float(self.num_clusters()) / 0.05

    def do_compute(self, iteration_result, ref_matrix):
        """compute method
        Note: will return None if not computed yet and the result of a previous
        scoring if the function is not supposed to actually run in this iteration
        """
        global SET_MATRIX, SET_MEMBERSHIP, SET_SET_TYPE
        logging.info("Compute scores for set enrichment...")
        start_time = util.current_millis()
        matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                               self.gene_names())
        use_multiprocessing = self.config_params[scoring.KEY_MULTIPROCESSING]
        SET_MATRIX = self.ratios
        SET_MEMBERSHIP = self.membership
        ref_min_score = ref_matrix.min()
        logging.info('REF_MIN_SCORE: %f', ref_min_score)

        set_filepath = os.path.join(self.config_params['output_dir'],
                                    'setEnrichment_set.csv')
        pval_filepath = os.path.join(self.config_params['output_dir'],
                                     'setEnrichment_pvalue.csv')

        for set_type in self.__set_types:
            SET_SET_TYPE = set_type
            logging.info("PROCESSING SET TYPE '%s'", set_type.name)
            start1 = util.current_millis()
            if use_multiprocessing:
                with util.get_mp_pool(self.config_params) as pool:
                    results = pool.map(compute_cluster_score,
                                       [(cluster, self.bonferroni_cutoff(), ref_min_score)
                                        for cluster in xrange(1, self.num_clusters() + 1)])
            else:
                results = []
                for cluster in xrange(1, self.num_clusters() + 1):
                    results.append(compute_cluster_score((cluster, self.bonferroni_cutoff(), ref_min_score)))

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
                    matrix.values[row][cluster - 1] = scores[row]
            setFile.write('\n'+str(iteration_result['iteration'])+','+','.join([str(i) for i in minSets]))
            pvFile.write('\n'+str(iteration_result['iteration'])+','+','.join([str(i) for i in pValues]))
            setFile.close()
            pvFile.close()

        logging.info("SET ENRICHMENT FINISHED IN %f s.\n",
                     (util.current_millis() - start_time) / 1000.0)
        return matrix

    def run_logs(self):
        """return the run logs"""
        return [self.run_log]


def compute_cluster_score(args):
    """Computes the cluster score for a given set type"""
    cluster, cutoff, ref_min_score = args
    set_type = SET_SET_TYPE
    matrix = SET_MATRIX
    cluster_rows = sorted(SET_MEMBERSHIP.rows_for_cluster(cluster))
    cluster_genes = [gene for gene in cluster_rows if gene in set_type.genes()]
    overlap_sizes = []
    set_sizes = []

    for set_name in sorted(set_type.sets.keys()):
        eset = set_type.sets[set_name]
        set_genes = eset.genes_above_cutoff()
        set_sizes.append(len(set_genes))

        intersect = set(cluster_genes).intersection(set_genes)
        overlap_sizes.append(len(intersect))

    num_sets = len(set_type.sets)
    num_genes = len(set_type.genes())
    phyper_n = list(np.array([num_genes] * num_sets) - np.array(set_sizes))
    phyper_k = [len(cluster_genes)] * num_sets

    enrichment_pvalues = np.array(util.phyper(overlap_sizes, set_sizes, phyper_n, phyper_k))
    min_pvalue = min(enrichment_pvalues[np.isfinite(enrichment_pvalues)])
    min_index = np.where(enrichment_pvalues == min_pvalue)[0][0]
    min_set = sorted(set_type.sets.keys())[min_index]
    min_set_overlap = overlap_sizes[min_index]

    scores = np.zeros(matrix.num_rows)
    if min_set_overlap > 0:
        min_genes = set_type.sets[min_set].genes()
        min_genes = [gene for gene in min_genes
                     if gene in matrix.row_names]
        min_indexes = matrix.row_indexes_for(min_genes)

        if set_type.sets[min_set].cutoff == 'discrete':
            overlap_genes = set(cluster_genes).intersection(set(min_genes))
            overlap_indexes = matrix.row_indexes_for(overlap_genes)
            scores[min_indexes] = 0.5
            scores[overlap_indexes] = 1.0
        else:
            # scaling
            for index in min_indexes:
                scores[index] = set_type.sets[min_set].weights[index]
            min_weight = min(min_set_weights)
            max_weight = max(min_set_weights)
            weight_range = max_weight - min_weight
            scores[min_indexes] -= min_weight
            scores[min_indexes] /= weight_range

        if min_pvalue <= cutoff:
            dampened_pvalue = 1
        else:
            dampened_pvalue = math.log10(min_pvalue) / math.log10(cutoff)

        scores[scores != 0.0] = dampened_pvalue / scores[scores != 0.0]
        scores *= ref_min_score

    return scores, min_set, min_pvalue
