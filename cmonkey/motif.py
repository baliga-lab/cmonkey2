"""motif.py - cMonkey motif related processing
This module captures the motif-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import logging
import numpy
import membership as memb
import datamatrix as dm


DISTANCE_UPSTREAM_SEARCH = (-20, 150)  # used to select sequences
DISTANCE_UPSTREAM_SCAN = (-30, 250)    # used for background distribution
MIN_CLUSTER_ROWS_ALLOWED = 3
MAX_CLUSTER_ROWS_ALLOWED = 70


# Applicable sequence filters
def unique_filter(seqs, feature_ids, distance):
    """returns a map that contains only the keys that are in
    feature_ids and only contains unique sequences"""
    unique_seqs = {}
    for feature_id in feature_ids:
        if (feature_id in seqs
            and seqs[feature_id] not in unique_seqs.values()):
            unique_seqs[feature_id] = seqs[feature_id]
    return unique_seqs


def get_remove_low_complexity_filter(meme_suite):
    """Factory method that returns a low complexity filter"""
    def remove_low_complexity(seqs, feature_ids, distance):
        """low-complexity filter that depends on meme"""
        return meme_suite.remove_low_complexity(seqs)
    return remove_low_complexity


def remove_atgs_filter(seqs, feature_ids, distance):
    """a filter removes the ATG's from the sequence, this
    just masks a window of 4 letters with N's"""
    for feature_id in seqs:
        chars = [c for c in seqs[feature_id]]
        chars[distance[1]:distance[1] + 4] = "NNNN"
        seqs[feature_id] = "".join(chars)
    return seqs


def compute_scores(meme_suite, organism, membership,
                   used_sequences,
                   distance, sequence_filters, pvalue_filter):
    """Compute motif scores. In order to influence the sequences
    that go into meme, the user can specify a list of sequence filter
    functions that have the signature
    (seqs, feature_ids, distance) -> seqs
    These filters are applied in the order they appear in the list
    """
    MIN_LOG_SCORE = -20.0

    def apply_sequence_filters(filters, seqs, feature_ids, distance):
        """apply all filters in the filters list in order"""
        for sequence_filter in filters:
            seqs = sequence_filter(seqs, feature_ids, distance)
        return seqs
    cluster_pvalues = {}

    for cluster in range(1, membership.num_clusters() + 1):
        logging.info("compute motif scores for cluster %d", cluster)
        genes = sorted(membership.rows_for_cluster(cluster))
        feature_ids = organism.feature_ids_for(genes)
        seqs = organism.sequences_for_genes(genes, distance, upstream=True)
        seqs = apply_sequence_filters(sequence_filters, seqs, feature_ids,
                                      distance)
        if (len(seqs) >= MIN_CLUSTER_ROWS_ALLOWED
            and len(seqs) <= MAX_CLUSTER_ROWS_ALLOWED):
            #logging.info("# seqs (= %d) within limits, continue " +
            #             "processing, seqs are: %s",
            #             len(seqs), str(seqs))
            pe_values, annotations = meme_suite.run_meme(seqs, used_sequences)

            pvalues = {}
            for feature_id, pvalue, evalue in pe_values:
                pvalues[feature_id] = numpy.log(pvalue)
            pvalues = pvalue_filter(pvalues)
            cluster_pvalues[cluster] = pvalues
            #for feature_id in pvalues:
            #    print "%s[%d] => %f" % (feature_id, cluster,
            #                            pvalues[feature_id])
            #print "PE-VALUES, CLUSTER: ", cluster
            #for feature_id, pvalue, evalue in pe_values:
            #    print "%s\t%f\t%f" % (feature_id, pvalue, evalue)

            #print "COUNT = %d" % len(pe_values)
            #print "------------------------"
            #print "ANNOTATIONS, CLUSTER: ", cluster
            #count = 0
            #for feature_id, annotation in annotations.items():
            #    for elem in annotation:
            #        print "%s\t%s" % (feature_id, str(elem))
            #        count += 1
            #print "COUNT = %d" % count
        else:
            logging.info("# seqs (= %d) outside of defined limits, skipping " +
                         "cluster %d", len(seqs), cluster)
        #print "CLUSTER PVALUES:"
        #print cluster_pvalues
    return cluster_pvalues


def make_min_value_filter(min_value):
    """A function generator which creates a value filter"""

    def min_value_filter(values):
        """all values that are below min_value are set to the
        minimum value above min_value"""
        allowed_vals = [value for key, value in values.items()
                        if value > min_value]
        result = {}
        if len(allowed_vals) == 0:
            logging.warn("all values are below the threshold -> no result !!")
        else:
            min_allowed = min(allowed_vals)
            for key, value in values.items():
                if value < min_value:
                    result[key] = min_allowed
                else:
                    result[key] = value
        return result
    return min_value_filter


class ScoringFunction(memb.ScoringFunctionBase):
    """Scoring function for motifs"""

    def __init__(self, organism, membership, matrix,
                 meme_suite, sequence_filters, pvalue_filter,
                 weight_func=None, interval=0):
        """creates a ScoringFunction"""
        memb.ScoringFunctionBase.__init__(self, membership,
                                          matrix, weight_func)
        self.__organism = organism
        self.__meme_suite = meme_suite
        self.__sequence_filters = sequence_filters
        self.__pvalue_filter = pvalue_filter

        # precompute the sequences for all genes that are referenced in the
        # input ratios, they are used as a basis to compute the background
        # distribution for every cluster
        self.__used_seqs = organism.sequences_for_genes(
            sorted(matrix.row_names()),
            DISTANCE_UPSTREAM_SCAN,
            upstream=True)
        self.__interval = interval
        self.__reverse_map = self.__build_reverse_map(matrix)

    def __build_reverse_map(self, matrix):
        """build a map that reconstructs the original row name from
        a feature id"""
        def feature_id_for(gene):
            """convenience method to return the feature id for a gene"""
            feature_ids = self.__organism.feature_ids_for([gene])
            if len(feature_ids) > 0:
                return feature_ids[0]
            else:
                return None

        result = {}
        for row_name in matrix.row_names():
            feature_id = feature_id_for(row_name)
            if feature_id != None:
                result[feature_id] = row_name
        return result

    def compute(self, iteration):
        """compute method, iteration is the 0-based iteration number"""
        if (self.__interval == 0 or
            (iteration > 0 and (iteration % self.__interval == 0))):
            print "RUN MOTIF SCORING IN ITERATION ", iteration
            pvalues = compute_scores(self.__meme_suite,
                                     self.__organism,
                                     self.membership(),
                                     self.__used_seqs,
                                     DISTANCE_UPSTREAM_SEARCH,
                                     self.__sequence_filters,
                                     self.__pvalue_filter)
            remapped = {}
            for cluster in pvalues:
                pvalues_k = pvalues[cluster]
                pvalues_genes = {}
                for feature_id, pvalue in pvalues_k.items():
                    pvalues_genes[self.__reverse_map[feature_id]] = pvalue
                remapped[cluster] = pvalues_genes
            #print remapped
            # convert remapped to an actual scoring matrix
            matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                                   self.gene_names())
            for row in matrix.row_names():
                row_index = matrix.row_names().index(row)
                for cluster in range(1, self.num_clusters() + 1):
                    if (cluster in remapped.keys() and
                        row in remapped[cluster].keys()):
                        matrix[row_index][cluster - 1] = remapped[cluster][row]
            return matrix
        else:
            return None
