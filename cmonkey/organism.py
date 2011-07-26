"""organism.py - organism-specific functionality in cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
from util import DelimitedFile, DelimitedFileMapper, best_matching_links
import re
import logging
import thesaurus
from seqtools import extract_upstream


def make_kegg_code_mapper(dfile):
    """returns a function that maps an organism code to a KEGG organism
    name"""
    return DelimitedFileMapper(dfile, 1, 3).lookup


def make_go_taxonomy_mapper(dfile):
    """returns a function that maps an RSAT organism name to a GO
    taxonomy id"""
    return DelimitedFileMapper(dfile, 0, 1).lookup


class RsatSpeciesInfo:  # pylint: disable-msg=R0903
    """A class to store species information retrieved from an RSAT database
    mirror. This is a mere value object"""

    def __init__(self, rsatdb, species, is_eukaryote, taxonomy_id):
        """create an instance of RsatSpeciesInfo"""
        self.rsatdb = rsatdb
        self.species = species
        self.is_eukaryote = is_eukaryote
        self.taxonomy_id = taxonomy_id


class Feature:  # pylint: disable-msg=R0902,R0903
    """representation of a feature. Just a value object"""

    def __init__(self, feature_id, feature_type, name, contig,
                 start, end, reverse):
        """Create a Feature instance"""
        # pylint: disable-msg=R0913
        self.__feature_id = feature_id
        self.__feature_type = feature_type
        self.__name = name
        self.__contig = contig
        self.__start = start
        self.__end = end
        self.__reverse = reverse

    def contig(self):
        """returns this feature's contig name"""
        return self.__contig

    def start(self):
        """returns this feature's start position on the strand"""
        return self.__start

    def end(self):
        """returns this feature's end position on the strand"""
        return self.__end

    def is_reverse(self):
        """returns whether feature is on the reverse strand"""
        return self.__reverse


def make_rsat_organism_mapper(rsatdb):
    """return a function that maps from a KEGG organism name to
    related RSAT information
    """
    def is_eukaryote(rsat_organism):
        """determine whether this organism is an eukaryote"""
        organism_text = rsatdb.get_organism(rsat_organism)
        return re.search('Eukaryota', organism_text) != None

    def get_taxonomy_id(rsat_organism):
        """Determine the taxonomy data from the RSAT database"""
        organism_names_dfile = DelimitedFile.create_from_text(
            rsatdb.get_organism_names(rsat_organism), comment='--')
        return organism_names_dfile.lines()[0][0]

    def mapper_fun(kegg_organism):
        """Mapper function to return basic information about an organism
        stored in the RSAT database. Only the genes in gene_names will
        be considered in the construction"""
        rsat_organism = best_matching_links(
            kegg_organism,
            rsatdb.get_directory())[0].rstrip('/')
        return RsatSpeciesInfo(rsatdb, rsat_organism,
                               is_eukaryote(rsat_organism),
                               get_taxonomy_id(rsat_organism))
    return mapper_fun


class OrganismFactory:
    """Factory to create an organism. Construction of an organism
    instance is relatively complex and costly, so it is coordinated
    here. Information has to be pulled together from various databases
    which are provided to the factory as configuration parameters.
    Note: this factory is biased towards microbial organisms and
    pulls information from
    - RSAT
    - STRING
    - GO
    - Microbes Online
    For other types of organisms, a different factory should be used
    """

    def __init__(self, code2kegg_organism,
                 rsat_organism_info,
                 get_go_taxonomy_id,
                 network_factories):
        """create a OrganismFactory instance"""
        self.code2kegg_organism = code2kegg_organism
        self.rsat_organism_info = rsat_organism_info
        self.get_taxonomy_id = get_go_taxonomy_id
        self.__network_factories = network_factories

    def create(self, organism_code):
        """factory method to create an organism from a code"""
        logging.info("Creating organism object for code '%s'...",
                     organism_code)
        kegg_organism = self.code2kegg_organism(organism_code)
        logging.info('KEGG organism: %s', kegg_organism)
        rsat_info = self.rsat_organism_info(kegg_organism)
        logging.info('RSAT info retrieved: %s', rsat_info.species)
        go_taxonomy_id = self.get_taxonomy_id(
            rsat_info.species.replace('_', ' '))
        logging.info('GO taxonomy id: %s', str(go_taxonomy_id))
        return Organism(organism_code, kegg_organism, rsat_info,
                        go_taxonomy_id, self.__network_factories)


class Organism:
    """Abstraction of an organism in cMonkey. It captures all organism-specific
    aspects. For now, we assume microbes only, but keep the interface generic
    so the algorithm will work on any type of organism"""

    def __init__(self, code, kegg_organism, rsat_info,
                 go_taxonomy_id, network_factories):
        """create an Organism instance"""
        self.code = code
        self.kegg_organism = kegg_organism
        self.__network_factories = network_factories
        self.__rsat_info = rsat_info
        self.go_taxonomy_id = go_taxonomy_id
        self.__synonyms = None  # lazy load
        self.__features = None
        self.__sequences = {}

    def species(self):
        """Retrieves the species of this object"""
        return self.__rsat_info.species

    def taxonomy_id(self):
        """Returns the taxonomy id"""
        return self.__rsat_info.taxonomy_id

    def is_eukaryote(self):
        """Determines whether this object is an eukaryote"""
        return self.__rsat_info.is_eukaryote

    def cog_organism(self):
        """returns the COG organism name"""
        return self.code.capitalize()

    def feature(self, gene_alias):
        """Returns this object's features"""
        return self.__features[self.feature_id_for(gene_alias)]

    def features_for_genes(self, gene_aliases):
        """returns a map of features for the specified list of genes aliases"""
        return self.__read_features_and_contigs(gene_aliases)[0]

    def init_with(self, gene_aliases, distance=(-30, 250)):
        """initialize this Organism's genome for usage in cMonkey"""
        self.__features, contigs = self.__read_features_and_contigs(
            gene_aliases)
        logging.info("Contigs: %s", str(contigs))
        logging.info("# Features read: %d", len(self.__features))
        self.__read_sequences(contigs, self.__features, distance)

    def feature_id_for(self, gene_alias):
        """retrieve the original feature id for a gene alias"""
        return self.__thesaurus()[gene_alias]

    def __thesaurus(self):
        """reads the thesaurus from a feature_names file"""
        if not self.__synonyms:
            feature_names_dfile = DelimitedFile.create_from_text(
                self.__rsatdb().get_feature_names(self.species()),
                comment='--')
            self.__synonyms = thesaurus.create_from_rsat_feature_names(
                feature_names_dfile, [thesaurus.strip_vng_modification])
        return self.__synonyms

    def __read_features_and_contigs(self, gene_names):
        """Reads RSAT features from a feature.tab file and returns a
        dictionary with feature ids as keys only the features that
        are in gene_names are actually read"""

        def add_feature_and_contig(features, contigs, feature_id, line):
            """Creates and adds a feature and associated contig from current
            DelimitedFile line"""
            contig = line[3]
            is_reverse = False
            if line[6] == 'R':
                is_reverse = True

            features[feature_id] = Feature(feature_id, line[1],
                                           line[2],
                                           contig,
                                           int(line[4]), int(line[5]),
                                           is_reverse)
            if contig not in contigs:
                contigs.append(contig)

        features = {}
        contigs = []
        synonyms = self.__thesaurus()
        id_names = [synonyms[name] for name in gene_names if name in synonyms]

        dfile = DelimitedFile.create_from_text(
            self.__rsatdb().get_features(self.species()), comment='--')
        for line in dfile.lines():
            feature_id = line[0]
            if feature_id in id_names:
                add_feature_and_contig(features, contigs, feature_id, line)

        return (features, contigs)

    def __rsatdb(self):
        """internal method to return the RSAT db link"""
        return self.__rsat_info.rsatdb

    def __read_sequences(self, contigs, features, distance):
        """for each feature, extract and set its sequence"""
        contig_seqs = {}
        for contig in contigs:
            contig_seqs[contig] = self.__rsatdb().get_contig_sequence(
                self.species(), contig)

        for feature_id in features:
            feature = features[feature_id]
            self.__sequences[feature_id] = extract_upstream(
                contig_seqs[feature.contig()],
                feature.start(),
                feature.end(),
                feature.is_reverse(),
                distance)

    def networks(self):
        """return the networks that can be generated by this
        organism"""
        result = []
        for make_network in self.__network_factories:
            result.append(make_network(self))
        return result

    def __str__(self):
        result = "Organism Type: %s\n" % self.__class__.__name__
        result += (("Code: '%s'\nKEGG: '%s'\nRSAT: '%s'\nCOG: '%s'\n" +
                   "GO Taxonomy Id: %s\n") %
                   (self.code, self.kegg_organism, self.__rsat_info.species,
                    self.cog_organism(), self.go_taxonomy_id))
        return result


__all__ = ['make_kegg_code_mapper', 'make_go_taxonomy_mapper',
           'make_rsat_organism_mapper', 'subsequence'
           'Organism', 'OrganismFactory']
