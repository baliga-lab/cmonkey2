"""set_enrichment.py - cMonkey set_enrichment scoring.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import util
import math
import logging
import scoring
import numpy as np
import rpy2.robjects as robjects
import datamatrix as dm


class EnrichmentSet:
    """Enrichment set representation"""
    def __init__(self, cutoff):
        """instance creation"""
        self.genes = []
        self.weights = []
        self.cutoff = cutoff
        self.__genes_above_cutoff = None

    def add(self, elem, weight):
        """Adds a gene and a weight to this set"""
        self.genes.append(elem)
        self.weights.append(weight)

    def genes_above_cutoff(self):
        """returns the genes that have a weight above the cutoff"""
        if self.__genes_above_cutoff == None:
            self.__genes_above_cutoff = []
            for index in xrange(0, len(self.genes)):
                if self.cutoff == 'discrete':
                    self.__genes_above_cutoff.append(self.genes[index])
                elif self.weights[index] >= self.cutoff:
                    self.__genes_above_cutoff.append(self.genes[index])
        return self.__genes_above_cutoff

    def __repr__(self):
        return ("Enrichment set: Cutoff = %s, # genes: %d # above cutoff: %d" %
                (self.cutoff, len(self.genes), len(self.genes_above_cutoff())))


class SetType:
    """Set type representation"""

    def __init__(self, name, sets):
        """instance creation"""
        self.name = name
        self.sets = sets
        self.__genes = None

    def genes(self):
        """All genes contained in the sets"""
        if self.__genes == None:
            self.__genes = set()
            for enrichment_set in self.sets.values():
                for gene in enrichment_set.genes:
                    self.__genes.add(gene)
        return self.__genes

    def __repr__(self):
        """string representation"""
        result = "SetType['%s'] = {\n" % self.name
        for key, value in self.sets.items():
            result += "%s -> %s\n" % (key, value)
        result += "}"
        return result

    @classmethod
    def read_csv(cls, name, infile, cutoff=None, sep=','):
        """reads a set from a CSV file"""
        dfile = util.DelimitedFile.read(infile, sep)
        sets = {}
        for line in dfile.lines():
            if line[0] not in sets:
                sets[line[0]] = EnrichmentSet('discrete')
            sets[line[0]].add(line[1].upper(), 1)
        return SetType(name, sets)


class ScoringFunction(scoring.ScoringFunctionBase):
    """Network scoring function"""

    def __init__(self, membership, matrix, set_types,
                 weight_func=None,
                 interval=0, config_params=None):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, weight_func,
                                             config_params)
        self.__interval = interval
        self.__set_types = set_types
        # stores (min_set, pvalue) pairs for each cluster and set type
        # for the last run of the function
        self.__last_min_enriched_set = {}
        for set_type in set_types:
            self.__last_min_enriched_set[set_type] = {}

    def bonferroni_cutoff(self):
        """Bonferroni cutoff value"""
        return float(self.num_clusters()) / 0.05

    def compute(self, iteration, ref_matrix):
        """compute method"""
        if (self.__interval == 0 or
            (iteration > 0 and (iteration % self.__interval == 0))):
            logging.info("Compute scores for set enrichment...")
            start_time = util.current_millis()
            matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                                   self.gene_names())
            for set_type in self.__set_types:
                #logging.info("PROCESSING SET TYPE [%s]", repr(set_type))
                for cluster in xrange(1, self.num_clusters() + 1):
                    scores = self.__compute_cluster_score(set_type,
                                                          cluster,
                                                          ref_matrix)
                    for row in xrange(len(self.gene_names())):
                        matrix[row][cluster - 1] = scores[row]
            logging.info("SET ENRICHMENT FINISHED IN %f s.\n",
                         (util.current_millis() - start_time) / 1000.0)
            return matrix
        else:
            return None

    def __compute_cluster_score(self, set_type, cluster, ref_matrix):
        """Computes the cluster score for a given set type"""
        cluster_rows = self.rows_for_cluster(cluster)
        cluster_genes = [gene for gene in cluster_rows
                         if gene in set_type.genes()]
        overlap_sizes = []
        set_sizes = []

        for set_name, eset in set_type.sets.items():
            set_genes = eset.genes_above_cutoff()
            intersect = set(cluster_genes).intersection(set_genes)
            overlap_sizes.append(len(intersect))
            set_sizes.append(len(set_genes))

        num_sets = len(set_type.sets)
        phyper_n = (np.array([len(set_type.genes())
                              for _ in xrange(num_sets)]) -
                    np.array(set_sizes))
        phyper_n = [value for value in phyper_n]
        phyper_k = [len(cluster_genes) for _ in xrange(num_sets)]
        enrichment_pvalues = list(util.phyper(overlap_sizes, set_sizes,
                                              phyper_n, phyper_k))
        min_pvalue = min(enrichment_pvalues)
        min_index = enrichment_pvalues.index(min_pvalue)
        min_set = set_type.sets.keys()[min_index]
        min_set_overlap = overlap_sizes[min_index]
        if min_set_overlap > 0:
            scores = [0.0 for _ in xrange(self.matrix().num_rows())]
            min_genes = set_type.sets[min_set].genes
            min_genes = [gene for gene in min_genes
                         if gene in self.gene_names()]
            min_indexes = self.matrix().row_indexes(min_genes)

            if set_type.sets[min_set].cutoff == 'discrete':
                overlap_genes = set(cluster_genes).intersection(set(min_genes))
                overlap_indexes = self.matrix().row_indexes(overlap_genes)
                for index in min_indexes:
                    scores[index] = 0.5
                for index in overlap_indexes:
                    scores[index] = 1.0
            else:
                min_set_weights = []
                for index in min_indexes:
                    min_set_weights.append(
                        set_type.sets[min_set].weights[index])
                min_weight = min(min_set_weights)
                max_weight = max(min_set_weights)
                for index in min_indexes:
                    scores[index] = min_set_weights[index] - min_weight
                    scores[index] = min_set_weights[index] / max_weight

            dampened_pvalue = enrichment_pvalues[min_index]
            if dampened_pvalue <= self.bonferroni_cutoff():
                dampened_pvalue = 1
            else:
                dampened_pvalue = (math.log10(dampened_pvalue) /
                                   math.log10(self.bonferroni_cutoff()))
            scores = [dampened_pvalue / score if score != 0.0 else score
                      for score in scores]
            min_ref_score = ref_matrix.min()
            scores = [score * min_ref_score for score in scores]
        else:
            scores = [0.0 for _ in xrange(self.matrix().num_rows())]

        # store the best enriched set determined
        self.__last_min_enriched_set[set_type][cluster] = (min_set, min_pvalue)
        return scores
