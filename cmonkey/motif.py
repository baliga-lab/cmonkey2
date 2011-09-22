"""motif.py - cMonkey motif related processing
This module captures the motif-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import logging
import numpy


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
    MIN_LOG_SCORE = -20.0

    """Compute motif scores. In order to influence the sequences
    that go into meme, the user can specify a list of sequence filter
    functions that have the signature
    (seqs, feature_ids, distance) -> seqs
    These filters are applied in the order they appear in the list
    """
    def apply_sequence_filters(filters, seqs, feature_ids, distance):
        """apply all filters in the filters list in order"""
        for sequence_filter in filters:
            seqs = sequence_filter(seqs, feature_ids, distance)
        return seqs

    for cluster in range(1, membership.num_clusters() + 1):
        logging.info("compute motif scores for cluster %d", cluster)
        genes = sorted(membership.rows_for_cluster(cluster))
        feature_ids = organism.feature_ids_for(genes)
        seqs = organism.sequences_for_genes(genes, distance, upstream=True)
        seqs = apply_sequence_filters(sequence_filters, seqs, feature_ids,
                                      distance)
        if (len(seqs) >= MIN_CLUSTER_ROWS_ALLOWED
            and len(seqs) <= MAX_CLUSTER_ROWS_ALLOWED):
            logging.info("# seqs (= %d) within limits, continue processing, " +
                         "seqs are: %s",
                         len(seqs), str(seqs))
            pe_values, annotations = meme_suite.run_meme(seqs, used_sequences)

            pvalues = {}
            for feature_id, pvalue, evalue in pe_values:
                pvalues[feature_id] = numpy.log(pvalue)
            pvalues = pvalue_filter(pvalues)
            for feature_id in pvalues:
                print "%s[%d] => %f" % (feature_id, cluster,
                                        pvalues[feature_id])
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


class ScoringFunction:
    """Scoring function for motifs"""

    def __init__(self, organism, membership, meme_suite, distance,
                 used_seqs, sequence_filters, pvalue_filter):
        """creates a ScoringFunction"""
        self.__organism = organism
        self.__membership = membership
        self.__meme_suite = meme_suite
        self.__distance = distance
        self.__used_seqs = used_seqs  # TODO: make this field computed
        self.__sequence_filters = sequence_filters
        self.__pvalue_filter = pvalue_filter

    def compute(self, matrix):
        """compute method"""
        return compute_scores(self.__meme_suite,
                              self.__organism,
                              self.__membership,
                              self.__used_seqs,
                              self.__distance,
                              self.__sequence_filters,
                              self.__pvalue_filter)
