# vi: sw=4 ts=4 et:
"""organism.py - organism-specific functionality in cMonkey
This module captures a microbial organism that receives data
from Microbes Online and RSAT

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import string
import logging
import thesaurus
import util
import seqtools as st
import microbes_online as mo
import collections
import patches


def make_kegg_code_mapper(dfile):
    """returns a function that maps an organism code to a KEGG organism
    name"""
    return util.DelimitedFileMapper(dfile, 1, 3).__getitem__


def make_go_taxonomy_mapper(dfile):
    """returns a function that maps an RSAT organism name to a GO
    taxonomy id"""
    return util.DelimitedFileMapper(dfile, 0, 1).__getitem__


class RsatSpeciesInfo:
    """RSAT description of the organism"""
    def __init__(self, rsatdb, kegg_organism, species, taxonomy_id):
        """determine RSAT information using the RSAT database object"""
        self.__rsatdb = rsatdb

        # in many cases, the fuzzy match delivers the correct RSAT organism
        # name, but there are exceptions
        if kegg_organism in patches.KEGG_EXCEPTIONS:
            kegg_organism = patches.KEGG_EXCEPTIONS[kegg_organism]

        if species is None:
            self.species = rsatdb.get_rsat_organism(kegg_organism)
        else:
            self.species = species

        logging.info("KEGG = '%s' -> RSAT = '%s'", kegg_organism, self.species)

        if taxonomy_id is None:
            self.taxonomy_id = rsatdb.get_taxonomy_id(self.species)
        else:
            self.taxonomy_id = taxonomy_id

    def get_features(self):
        return self.__rsatdb.get_features(self.species)

    def get_feature_names(self):
        return self.__rsatdb.get_feature_names(self.species)

    def get_contig_sequence(self, contig):
        return self.__rsatdb.get_contig_sequence(self.species, contig)

    def go_species(self):
        return self.species.replace('_', ' ')


KEGGExceptions = {'Pseudomonas aeruginosa PAO1': 'Pseudomonas aeruginosa',
                  'Campylobacter jejuni NCTC11168': 'Campylobacter jejuni'}


"""
def make_rsat_organism_mapper(rsatdb):
    def mapper_fun(kegg_organism, rsat_organism, ncbi_code=None):
        return RsatSpeciesInfo(rsatdb, kegg_organism, rsat_organism, ncbi_code)
    return mapper_fun
"""

class OrganismBase:
    """The organism base class contains functionality that is likely to
    be the same among Organism implementations"""

    def __init__(self, code, network_factories, ratios=None):
        """Initialize the base class instance"""
        self.code = code
        logging.info("Creating networks...")
        self.__networks = []
        for make_network in network_factories:
            self.__networks.append(make_network(self, ratios))
        logging.info("Finished creating networks.")

    def networks(self):
        """returns this organism's networks"""
        return self.__networks

    def thesaurus(self):
        """Returns a map containing the alias -> gene mappings"""
        raise Exception("please implement me")

    def feature_ids_for(self, gene_aliases):
        """Helper method to retrieve a list of feature_ids for the
        specified alias list"""
        synonyms = self.thesaurus()
        return [synonyms[alias] for alias in gene_aliases if alias in synonyms]


class DummyOrganism(OrganismBase):
    """minimal organism class for gene expression-only runs"""

    def __init__(self):
        OrganismBase.__init__(self, 0, [])

    def thesaurus(self):
        return {}

    def species(self):
        return "Dummy organism"


class RSATOrganism(OrganismBase):
    """An organism class that does most things automatically by relying on information
    stored retrieved from RSAT."""

    # pylint: disable-msg=R0913,R0902
    def __init__(self, code, kegg_organism, rsat_info, go_taxonomy_id,
                 network_factories, search_distances, scan_distances,
                 use_operons=True, ratios=None):
        """create an Organism instance"""
        # microbe-specific network factories need access to synonyms
        # and rsat info, so initialize them here before the base class
        # init
        self.__synonyms = None  # lazy loaded
        self.__rsat_info = rsat_info
        self.use_operons = use_operons
        OrganismBase.__init__(self, code, network_factories, ratios=ratios)
        self.kegg_organism = kegg_organism
        self.go_taxonomy_id = go_taxonomy_id
        self.search_distances = search_distances
        self.scan_distances = scan_distances

    def species(self):
        """Retrieves the species of this object"""
        return self.__rsat_info.species

    def taxonomy_id(self):
        """Returns the taxonomy id"""
        return self.__rsat_info.taxonomy_id

    def cog_organism(self):
        """returns the COG organism name"""
        return self.code.capitalize()

    def thesaurus(self):
        """reads the thesaurus from a feature_names file. The thesaurus
        is also cached, because it is used many times
        """
        if not self.__synonyms:
            feature_names_dfile = util.dfile_from_text(
                self.__rsat_info.get_feature_names(),
                comment='--')
            self.__synonyms = thesaurus.create_from_rsat_feature_names(
                feature_names_dfile, [thesaurus.strip_vng_modification])
        return self.__synonyms

    def features_for_genes(self, genes):
        """returns a map of features for the specified list of genes aliases
        used for operon information"""
        return util.ThesaurusBasedMap(
            self.thesaurus(),
            self.read_features(self.feature_ids_for(genes)))

    def read_features(self, feature_ids):
        """Returns a list containing the features for the specified feature
        ids"""

        def read_feature(line):
            """Creates and adds a feature and associated contig from current
            DelimitedFile line"""
            contig = line[3]
            is_reverse = False
            if line[6] == 'R':
                is_reverse = True

            # note that feature positions can sometimes start with a '>'
            # or '<', so make sure it is stripped away
            return st.Feature(line[0], line[1], line[2],
                              st.Location(contig,
                                          int(string.lstrip(line[4], '<>')),
                                          int(string.lstrip(line[5], '<>')),
                                          is_reverse))

        features = {}
        dfile = util.dfile_from_text(self.__rsat_info.get_features(), comment='--')
        for line in dfile.lines:
            feature_id = line[0]
            if feature_id in feature_ids:
                features[feature_id] = read_feature(line)
        return features

    def read_sequences(self, features, distance, extractor):
        """for each feature, extract and set its sequence"""
        def unique_contigs():
            """extract the unique contigs from the input features"""
            result = []
            for feature in features.values():
                if feature.location.contig not in result:
                    result.append(feature.location.contig)
            return result

        sequences = {}
        contig_seqs = {contig: self.__rsat_info.get_contig_sequence(contig)
                       for contig in unique_contigs()}

        for key, feature in features.items():
            location = feature.location
            sequences[key] = extractor(
                contig_seqs[location.contig], location, distance)
        if len(sequences) == 0:
            logging.error('No sequences read for %s!' % self.code)
        return sequences

    def __str__(self):
        result = "Organism Type: %s\n" % self.__class__.__name__
        result += (("Code: '%s'\nKEGG: '%s'\nRSAT: '%s'\nCOG: '%s'\n" +
                   "GO Taxonomy Id: %s\n") %
                   (self.code, self.kegg_organism, self.species(),
                    self.cog_organism(), self.go_taxonomy_id))
        return result


class Microbe(RSATOrganism):
    """The standard organism in cmonkey is Microbe. It builds on the
    data dependencies established in RSATOrganism and adds a Microbes Online
    dependency to retrieve possible operon information """

    # pylint: disable-msg=R0913,R0902
    def __init__(self, code, kegg_organism, rsat_info,
                 go_taxonomy_id, microbes_online_db,
                 network_factories,
                 search_distances, scan_distances,
                 use_operons=True, ratios=None):
        """create an Organism instance"""
        RSATOrganism.__init__(self, code, kegg_organism,
                              rsat_info, go_taxonomy_id, network_factories,
                              search_distances, scan_distances, use_operons, ratios)
        self.__microbes_online_db = microbes_online_db
        self.__operon_mappings = None  # lazy loaded

    def sequences_for_genes_search(self, genes, seqtype='upstream'):
        """The default sequence retrieval for microbes is to
        fetch their operon sequences"""
        return self.operon_shifted_seqs_for(genes,
                                            self.search_distances[seqtype])

    def sequences_for_genes_scan(self, genes, seqtype='upstream'):
        """The default sequence retrieval for microbes is to
        fetch their operon sequences"""
        return self.operon_shifted_seqs_for(genes,
                                            self.scan_distances[seqtype])

    def operon_shifted_seqs_for(self, gene_aliases, distance):
        """returns a map of the gene_aliases to the feature-
        sequence tuple that they are actually mapped to.
        """
        def do_operon_shift():
            """Extract the (gene, head) pairs that are actually used"""
            operon_map = self.operon_map()
            synonyms = self.thesaurus()
            shifted_pairs = []
            aliases_not_found = []
            operons_not_found = []
            for alias in gene_aliases:
                if alias in synonyms:
                    gene = synonyms[alias]
                    if gene in operon_map:
                        shifted_pairs.append((gene, operon_map[gene]))
                    else:
                        operons_not_found.append(alias)
                        shifted_pairs.append((gene, gene))

                else:
                    aliases_not_found.append(alias)
            return shifted_pairs

        def unique_sequences(operon_pairs):
            """Returns the unique sequences for the specified operon pairs"""
            unique_feature_ids = []
            for _, head in operon_pairs:
                if head not in unique_feature_ids:
                    unique_feature_ids.append(head)
            features = self.read_features(unique_feature_ids)
            return self.read_sequences(features, distance,
                                         st.extract_upstream)

        if self.use_operons:
            shifted_pairs = do_operon_shift()
        else:
            # if operons should not be used, we simply map
            # the gene heads to themselves
            synonyms = self.thesaurus()
            valid_genes = [synonyms[alias]
                           for alias in gene_aliases if alias in synonyms]
            shifted_pairs = [(gene, gene) for gene in valid_genes]

        unique_seqs = unique_sequences(shifted_pairs)
        return {gene: unique_seqs[head] for gene, head in shifted_pairs}

    def operon_map(self):
        """Returns the operon map for this particular organism.
        Microbes Online works on VNG names, but RSAT is working on
        feature ids, so this function also maps VNG names to feature ids"""
        if not self.__operon_mappings:
            pairs = mo.get_operon_pairs(self.__microbes_online_db, self)
            synonyms = self.thesaurus()
            self.__operon_mappings = {synonyms[gene]: synonyms[head] for head, gene in pairs}
        return self.__operon_mappings


__all__ = ['make_kegg_code_mapper', 'make_go_taxonomy_mapper',
           'Organism', 'OrganismFactory']
