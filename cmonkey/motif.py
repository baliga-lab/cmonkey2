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


def compute_scores(meme_suite, organism, membership, used_sequences):
    """compute motif scores"""
    def filter_sequences(meme_suite, sorted_feature_ids, seqs,
                         distance):
        """filter out redundant and low-complexity sequences"""
        unique_seqs = {}
        for feature_id in sorted_feature_ids:
            if (feature_id in seqs
                and seqs[feature_id] not in unique_seqs.values()):
                unique_seqs[feature_id] = seqs[feature_id]
        return remove_atgs(meme_suite.remove_low_complexity(unique_seqs),
                           distance)

    def remove_atgs(sequence_map, distance):
        for feature_id in sequence_map:
            chars = [c for c in sequence_map[feature_id]]
            chars[distance[1]:distance[1] + 4] = "NNNN"
            sequence_map[feature_id] = "".join(chars)
        return sequence_map

    distance = DISTANCE_UPSTREAM_SEARCH
    for cluster in [2]:
        genes = sorted(membership.rows_for_cluster(cluster))
        feature_ids = organism.feature_ids_for(genes)
        seqs = organism.sequences_for_genes(genes, distance, upstream=True)
        seqs = filter_sequences(meme_suite, feature_ids, seqs, distance)
        if (len(seqs) >= MIN_CLUSTER_ROWS_ALLOWED
            and len(seqs) <= MAX_CLUSTER_ROWS_ALLOWED):
            logging.info("# seqs (= %d) within limits, continue processing, seqs are: %s",
                         len(seqs), str(seqs))
            meme_suite.run_meme(seqs, used_sequences)
        else:
            logging.info("# seqs (= %d) outside of defined limits, skipping " +
                         "cluster %d", len(seqs), cluster)
