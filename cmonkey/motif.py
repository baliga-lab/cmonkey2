"""motif.py - cMonkey motif related processing
This module captures the motif-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


DISTANCE_UPSTREAM_SEARCH = (-20, 150)
DISTANCE_UPSTREAM_SCAN = (-30, 150)


"""
def __get_operons(organism, genes):
    operon_map = organism.operon_map()
    #print operon_map
    operon_pairs = []
    for gene in genes:
        if gene in operon_map:
            operon_pairs.append((gene, operon_map[gene]))
    print operon_pairs
    operons = []
    for gene, operon in operon_pairs:
        print "Searching Operon: %s" % operon
        features = organism.features_for_genes([operon])
        operons.append(features[features.keys()[0]])
    print operons
    opseqs = organism.sequences_for_features(operons)
    print opseqs.keys()
"""


def compute_scores(organism, membership):
    """compute motif scores"""
    genes = sorted(membership.rows_for_cluster(1))
    seqs = organism.sequences_for_genes(genes, DISTANCE_UPSTREAM_SEARCH)
