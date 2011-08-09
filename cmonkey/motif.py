"""motif.py - cMonkey motif related processing
This module captures the motif-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
DISTANCE_UPSTREAM_SEARCH = (-20, 150)
DISTANCE_UPSTREAM_SCAN = (-30, 150)


def compute_scores(organism, membership):
    """compute motif scores"""
    genes = sorted(membership.rows_for_cluster(1))
    print genes
    seqs = organism.sequences_for_genes(genes, DISTANCE_UPSTREAM_SEARCH)
    print seqs
