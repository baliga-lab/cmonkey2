"""motif.py - cMonkey motif related processing
This module captures the motif-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import logging


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


def compute_scores(meme_suite, organism, membership, used_sequences,
                   distance, sequence_filters):
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

    for cluster in [2]:
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
            meme_suite.run_meme(seqs, used_sequences)
        else:
            logging.info("# seqs (= %d) outside of defined limits, skipping " +
                         "cluster %d", len(seqs), cluster)
