# vi: sw=4 ts=4 et:
"""microbes_online.py - operon prediction module cMonkey
For organisms that have operons, cMonkey can read use operon predictions
to improve its scoring. This module defines access to prediction data through
various services.
Note: This module uses VNG names, mostly to be easily comparable
----- with the reference cMonkey implementation.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import logging

import cmonkey.util as util
import cmonkey.network as network
import cmonkey.patches as patches

MICROBES_ONLINE_BASE_URL = 'http://www.microbesonline.org'
MYSQL_HOST = 'pub.microbesonline.org'
MYSQL_USER = 'guest'
MYSQL_PASSWD = 'guest'
MYSQL_DB = 'genomics'


class MicrobesOnlineOperonFile:
    """access to Microbes online operon prediction by providing a file"""

    def __init__(self, path):
        self.path = path

    def get_operon_predictions_for(self, organism_id):
        with open(self.path) as infile:
            return infile.read()


class MicrobesOnline:
    """Interface to Microbes Online web service"""

    def __init__(self, cache_dir,
                 base_url=MICROBES_ONLINE_BASE_URL):
        """creates a MicrobesOnline service instance"""
        self.base_url = base_url
        self.cache_dir = cache_dir

    def get_operon_predictions_for(self, organism_id):
        """Retrieve operon predictions for the specified organism"""
        logging.info("MicrobesOnline.get_operon_predictions_for(%s)",
                     organism_id)
        url = '/'.join([self.base_url, 'operons',
                       'gnc%s.named' % str(organism_id)])
        cache_file = '/'.join([self.cache_dir,
                              'gnc%s.named' % str(organism_id)])
        return util.read_url_cached(url, cache_file).decode('utf-8')


# Operon calculation functionality.
# Conceptually, this could go into an Organism class that has operons.
# It is still placed within the microbes_online module because
#
# 1. it is strongly related to the structure of Microbes Online Data
# 2. computation of operon maps and networks is pretty complicated
#    and it benefits from separate testing

def make_operon_pairs(operon, features):
    """take an operon as a list of gene names, determines the head out of
    these gene names and generates a (head, gene) for each gene in the
    operon.
    features is a map from a gene alias/feature id to a Feature object
    The head is is determined as follows:
    1. retrieve the gene coordinates for each gene in the operon
    2. if most genes are on the forward strand, the head is the one
       with the lowest start position
    3. if most genes are on the reverse strand, the head is the one
       with the highest end position
    This function returns an empty result if
    1. the same amount of genes are on the forward and reverse strand
    2. the gene coordinates can't be retrieved
    """
    def get_reverse_head(feature_map):
        """determine reverse head of the operon"""
        max_gene = None
        max_end = 0
        for (gene, feature) in feature_map.items():
            if feature.location.end > max_end:
                max_end = feature.location.end
                max_gene = gene
        return max_gene

    def get_forward_head(feature_map):
        """determine forward head of the operon"""
        min_gene = None
        min_start = sys.maxsize
        for (gene, feature) in feature_map.items():
            if feature.location.start < min_start:
                min_start = feature.location.start
                min_gene = gene
        return min_gene

    feature_map = {}  # mapping from VNG name to feature
    num_reverse = 0

    # make sure we only take the genes that we have genomic information
    # for, and ignore the rest
    available_operon_genes = []
    for gene in operon:
        if gene in features.keys():
            available_operon_genes.append(gene)
        else:
            logging.warn("Microbes Online operon gene '%s' not found in " +
                         "RSAT features", gene)

    for gene in available_operon_genes:
        feature_map[gene] = features[gene]
        if feature_map[gene].location.reverse:
            num_reverse += 1

    num_total = len(available_operon_genes)
    if num_total > 0:
        percent_reverse = float(num_reverse) / float(num_total)
        if percent_reverse > 0.6:
            head = get_reverse_head(feature_map)
        elif percent_reverse < 0.4:
            head = get_forward_head(feature_map)
        else:
            logging.warning("can't determine head of operon - amounts " +
                            "of reverse and forward genes are too similar " +
                            "(%f-%f)",
                            percent_reverse, 1.0 - percent_reverse)
            return []
        return [(head, gene) for gene in available_operon_genes]
    else:
        logging.warning("Operon did not contain any available genes")
        return []


def build_operons(names1, names2):
    """build the list of operons given two name lists"""
    def first_row_containing(alist, aname):
        """returns in a list containing the specified name"""
        for row in alist:
            if aname in row:
                return row
        return None

    operons = []
    for i in range(len(names1)):
        found = first_row_containing(operons, names1[i])
        if found:
            found.append(names2[i])
        else:
            operons.append([names1[i], names2[i]])
    return operons


def __make_operons_from_predictions(predictions, organism):
    """returns a operon list and feature list for the
    specified predictions and organism"""

    def build_names():
        """builds the gene name lists from the predictions"""
        names1 = []
        names2 = []
        for prediction in predictions:
            if prediction[0] not in names1:
                names1.append(prediction[0])
            if prediction[1] not in names2:
                names2.append(prediction[1])
        return names1, names2

    names1, names2 = build_names()
    features = organism.features_for_genes(names1 + names2)
    operons = build_operons(names1, names2)
    logging.info("%d operons created", len(operons))
    return operons, features


def make_pairs_from_predictions(predictions, organism):
    """return a list of predictions into a list of network edges"""
    operons, features = __make_operons_from_predictions(predictions, organism)
    pairs = []
    for operon in operons:
        pairs.extend(make_operon_pairs(operon, features))
    return pairs


def get_operon_pairs(microbes_online, organism):
    """returns a list of (head, gene) pairs that were derived from
    an operon prediction file for an organism from Microbes Online
    Used for retrieving genes that have an operon shift
    """
    preds = __get_predictions(microbes_online, organism)
    return make_pairs_from_predictions(preds, organism)


def __get_predictions(microbes_online, organism):
    """reads the operon predictions for a given organism from MicrobesOnline"""
    preds_text = microbes_online.get_operon_predictions_for(
        organism.taxonomy_id())
    dfile = util.dfile_from_text(preds_text, has_header=True)
    code = organism.code
    preds = [(patches.patch_mo_gene(code, line[2]),
              patches.patch_mo_gene(code, line[3]))
             for line in dfile.lines if line[6] == 'TRUE']
    logging.info("%d prediction pairs read", len(preds))
    return preds


def get_network_factory(microbes_online, max_operon_size, weight):
    """function to create a network factory method"""

    def get_operon_edges(microbes_online, organism):
        """gets network edges"""
        operons, _ = __make_operons_from_predictions(
            __get_predictions(microbes_online, organism), organism)
        edges = []
        for operon in operons:
            if len(operon) <= max_operon_size:
                combs = util.kcombinations(operon, 2)
                edges.extend([(comb[0], comb[1], 1000.0)
                              for comb in combs if comb[0] != comb[1]])
            else:
                logging.warn("dropped operon from network (max_operon_size " +
                             "exceeded): %s", str(operon))
        return edges

    def make_network(organism, ratios=None, check_size=True):
        """factory method to create a network from operon predictions"""

        logging.info("MicrobesOnline - make_network()")
        edges = get_operon_edges(microbes_online, organism)
        logging.info("%d edges computed", len(edges))
        return network.Network.create('operons', edges, weight, organism,
                                      ratios, check_size)

    return make_network


#######################################################################
#### DIRECT MYSQL ACCESS
#######################################################################

def mysql_connect():
    """create a database connection"""
    return mysql.connect(host=MYSQL_HOST, user=MYSQL_USER,
                         passwd=MYSQL_PASSWD, db=MYSQL_DB)


def mysql_synonyms(conn, taxonomy_id):
    """Retrieve synonyms through direct MySQL access"""
    cursor = conn.cursor()
    cursor.execute('select distinct Locus.locusId, Synonym.name, ' +
                   'Synonym.type from ' +
                   'Scaffold join Locus using (scaffoldId) join ' +
                   'Synonym using (locusId, version) ' +
                   'where Locus.priority=1 and Locus.type=1 and ' +
                   'Scaffold.isActive=1 and Synonym.type in (0, 1, 3) and ' +
                   'Scaffold.taxonomyId= %(taxid)s order by ' +
                   'Locus.locusId, Synonym.type',
                   {'taxid': taxonomy_id})
    row = cursor.fetchone()
    while row:
        print ("%d: %s (%s)" % (row[0], row[1], row[2]))
        row = cursor.fetchone()


__all__ = ['MicrobesOnline', 'get_network_factory', 'get_operon_pairs']
