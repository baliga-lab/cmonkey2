"""organism.py - organism-specific functionality in cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
from util import DelimitedFile, DelimitedFileMapper, best_matching_links
import re
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

    def __init__(self, species, is_eukaryote, taxonomy_id,
                 features, contigs):
        """create an instance of RsatSpeciesInfo"""
        # pylint: disable-msg=R0913
        self.species = species
        self.is_eukaryote = is_eukaryote
        self.taxonomy_id = taxonomy_id
        self.features = features
        self.contigs = contigs


class Feature:  # pylint: disable-msg=R0902,R0903
    """representation of a feature. Just a value object"""

    def __init__(self, feature_id, feature_type, name, contig,
                 start, end, reverse):
        """Create a Feature instance"""
        # pylint: disable-msg=R0913
        self.feature_id = feature_id
        self.feature_type = feature_type
        self.name = name
        self.contig = contig
        self.start = start
        self.end = end
        self.reverse = reverse
        self.sequence = None

    def set_sequence(self, seq):
        """Sets this Feature's sequence"""
        self.sequence = seq


def make_rsat_organism_mapper(rsatdb, distance=(-30, 250)):
    """return a function that maps from a KEGG organism name to
    related RSAT information
    TODO: we also need to retrieve operon sequences if available
    """
    # pylint: disable-msg=R0912

    def read_rsat_features_and_contigs(rsat_organism, gene_names, synonyms):
        """Reads RSAT features from a feature.tab file and returns a
        dictionary with feature ids as keys only the features that
        are in gene_names are actually read"""
        features = {}
        contigs = []
        id_names = [synonyms[name] for name in gene_names if name in synonyms]

        dfile = DelimitedFile.create_from_text(
            rsatdb.get_features(rsat_organism), comment='--')
        for line in dfile.lines():
            feature_id = line[0]
            if feature_id in id_names:
                add_feature_and_contig(features, contigs, feature_id, line)

        return (features, contigs)

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

    def add_seqs_to_features(rsat_organism, contigs, features):
        """for each feature, extract and set its sequence"""
        contig_seqs = {}
        for contig in contigs:
            contig_seqs[contig] = rsatdb.get_contig_sequence(rsat_organism,
                                                             contig)
        for feature_id in features:
            feature = features[feature_id]
            feature.set_sequence(extract_upstream(contig_seqs[feature.contig],
                                                  feature.start, feature.end,
                                                  feature.reverse,
                                                  distance))

    def strip_vng_modification(gene):
        """strips 'm' modifier off a VNG name"""
        if re.match('VNG\d{4}.m$', gene):
            return gene.rstrip('m')
        else:
            return gene

    def read_synonyms(rsat_organism):
        """reads the thesaurus from a feature_names file"""
        feature_names_dfile = DelimitedFile.create_from_text(
            rsatdb.get_feature_names(rsat_organism), comment='--')
        return thesaurus.create_from_rsat_feature_names(
            feature_names_dfile, [strip_vng_modification])

    def is_eukaryote(rsat_organism):
        """determine whether this organism is an eukaryote"""
        organism_text = rsatdb.get_organism(rsat_organism)
        return re.search('Eukaryota', organism_text) != None

    def get_taxonomy_id(rsat_organism):
        """Determine the taxonomy data from the RSAT database"""
        organism_names_dfile = DelimitedFile.create_from_text(
            rsatdb.get_organism_names(rsat_organism), comment='--')
        return organism_names_dfile.lines()[0][0]

    def mapper_fun(kegg_organism, gene_names):
        """Mapper function to return basic information about an organism
        stored in the RSAT database. Only the genes in gene_names will
        be considered in the construction"""
        rsat_organism = best_matching_links(
            kegg_organism,
            rsatdb.get_directory())[0].rstrip('/')
        synonyms = read_synonyms(rsat_organism)
        features, contigs = read_rsat_features_and_contigs(rsat_organism,
                                                           gene_names,
                                                           synonyms)
        add_seqs_to_features(rsat_organism, contigs, features)
        # just debugging
        """
        for name in sorted(gene_names):
            if name in synonyms:
                feature = features[synonyms[name]]
                print("'%s' -> '%s' [%d, %d] contig: %s, rev: %s = '%s'" %
                      (name, synonyms[name], feature.start, feature.end,
                       feature.contig, str(feature.reverse), feature.sequence))
                       """
        return RsatSpeciesInfo(rsat_organism, is_eukaryote(rsat_organism),
                               get_taxonomy_id(rsat_organism),
                               features, contigs)
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
                 get_go_taxonomy_id):
        """create a OrganismFactory instance"""
        self.code2kegg_organism = code2kegg_organism
        self.rsat_organism_info = rsat_organism_info
        self.get_taxonomy_id = get_go_taxonomy_id

    def create(self, organism_code, gene_names):
        """factory method to create an organism from a code"""
        kegg_organism = self.code2kegg_organism(organism_code)
        rsat_info = self.rsat_organism_info(kegg_organism, gene_names)
        go_taxonomy_id = self.get_taxonomy_id(
            rsat_info.species.replace('_', ' '))
        if rsat_info.is_eukaryote:
            return Eukaryote(organism_code, kegg_organism, rsat_info,
                             go_taxonomy_id)
        else:
            return Prokaryote(organism_code, kegg_organism, rsat_info,
                              go_taxonomy_id)


class Organism:
    """Abstraction of an organism in cMonkey. It captures all organism-specific
    aspects"""

    def __init__(self, code, kegg_organism, rsat_info, go_taxonomy_id):
        """create an Organism instance"""
        self.code = code
        self.kegg_organism = kegg_organism
        self.rsat_info = rsat_info
        self.go_taxonomy_id = go_taxonomy_id

    def cog_organism(self):
        """returns the COG organism name"""
        return self.code.capitalize()

    def __str__(self):
        result = "Organism Type: %s\n" % self.__class__.__name__
        result += (("Code: '%s'\nKEGG: '%s'\nRSAT: '%s'\nCOG: '%s'\n" +
                   "GO Taxonomy Id: %s\nContigs: %s\n") %
                   (self.code, self.kegg_organism, self.rsat_info.species,
                    self.cog_organism(), self.go_taxonomy_id,
                    str(self.rsat_info.contigs)))
        return result


class Eukaryote(Organism):
    """Class to represent eukaryotes"""

    def is_eukaryote(self):  # pylint: disable-msg=R0201
        """always returns true"""
        return True


class Prokaryote(Organism):
    """Class to represent prokaryotes"""

    def is_eukaryote(self):  # pylint: disable-msg=R0201
        """always returns false"""
        return False


__all__ = ['make_kegg_code_mapper', 'make_go_taxonomy_mapper',
           'make_rsat_organism_mapper', 'subsequence'
           'Organism', 'OrganismFactory']
