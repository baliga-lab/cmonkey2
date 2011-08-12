"""motif.py - cMonkey motif related processing
This module captures the motif-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import logging


DISTANCE_UPSTREAM_SEARCH = (-20, 150)
DISTANCE_UPSTREAM_SCAN = (-30, 150)
MIN_CLUSTER_ROWS_ALLOWED = 3
MAX_CLUSTER_ROWS_ALLOWED = 70


def compute_scores(meme_suite, organism, membership):
    """compute motif scores"""
    def filter_sequences(meme_suite, sorted_feature_ids, seqs):
        """filter out redundant and low-complexity sequences"""
        unique_seqs = {}
        for feature_id in sorted_feature_ids:
            if (feature_id in seqs
                and seqs[feature_id] not in unique_seqs.values()):
                unique_seqs[feature_id] = seqs[feature_id]
        return meme_suite.remove_low_complexity(unique_seqs)

    for cluster in [2]:
        genes = sorted(membership.rows_for_cluster(cluster))
        feature_ids = organism.feature_ids_for(genes)
        seqs = organism.sequences_for_genes(genes, DISTANCE_UPSTREAM_SEARCH,
                                            upstream=True, motif_finding=True)
        seqs = filter_sequences(meme_suite, feature_ids, seqs)
        if (len(seqs) >= MIN_CLUSTER_ROWS_ALLOWED
            and len(seqs) <= MIN_CLUSTER_ROWS_ALLOWED):
            logging.info("# seqs (= %d) within limits, continue processing",
                         len(seqs))
        else:
            logging.info("# seqs (= %d) outside of defined limits, skipping " +
                         "cluster %d", len(seqs), cluster)
