# vi: sw=4 ts=4 et:
"""organism.py - organism-specific functionality in cMonkey
This module captures a microbial organism that receives data
from Microbes Online and RSAT

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import re
import string
import logging
import thesaurus
import util
import seqtools as st
import microbes_online as mo


def make_kegg_code_mapper(dfile):
    """returns a function that maps an organism code to a KEGG organism
    name"""
    return util.DelimitedFileMapper(dfile, 1, 3).__getitem__


def make_go_taxonomy_mapper(dfile):
    """returns a function that maps an RSAT organism name to a GO
    taxonomy id"""
    return util.DelimitedFileMapper(dfile, 0, 1).__getitem__


class RsatSpeciesInfo:  # pylint: disable-msg=R0903
    """A class to store species information retrieved from an RSAT database
    mirror. This is a mere value object"""

    def __init__(self, rsatdb, species, is_eukaryote, taxonomy_id):
        """create an instance of RsatSpeciesInfo"""
        self.rsatdb = rsatdb
        self.species = species
        self.is_eukaryote = is_eukaryote
        self.taxonomy_id = taxonomy_id


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
        organism_names_dfile = util.DelimitedFile.create_from_text(
            rsatdb.get_organism_names(rsat_organism), comment='--')
        return organism_names_dfile.lines()[0][0]

    def mapper_fun(kegg_organism):
        """Mapper function to return basic information about an organism
        stored in the RSAT database. Only the genes in gene_names will
        be considered in the construction"""
        rsat_organism = util.best_matching_links(
            kegg_organism,
            rsatdb.get_directory())[0].rstrip('/')
        return RsatSpeciesInfo(rsatdb, rsat_organism,
                               is_eukaryote(rsat_organism),
                               get_taxonomy_id(rsat_organism))
    return mapper_fun


class MicrobeFactory:
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

    # pylint: disable-msg=R0913
    def __init__(self, code2kegg_organism,
                 rsat_organism_info,
                 get_go_taxonomy_id,
                 microbes_online_db,
                 network_factories):
        """create a OrganismFactory instance"""
        self.__code2kegg_organism = code2kegg_organism
        self.__rsat_organism_info = rsat_organism_info
        self.__get_taxonomy_id = get_go_taxonomy_id
        self.__microbes_online_db = microbes_online_db
        self.__network_factories = network_factories

    def create(self, organism_code, search_distances,
               scan_distances):
        """factory method to create an organism from a code"""
        logging.info("Creating organism object for code '%s'...",
                     organism_code)
        kegg_organism = self.__code2kegg_organism(organism_code)
        logging.info('KEGG organism: %s', kegg_organism)
        rsat_info = self.__rsat_organism_info(kegg_organism)
        logging.info('RSAT info retrieved: %s', rsat_info.species)
        go_taxonomy_id = self.__get_taxonomy_id(
            rsat_info.species.replace('_', ' '))
        logging.info('GO taxonomy id: %s', str(go_taxonomy_id))
        return Microbe(organism_code, kegg_organism, rsat_info,
                       go_taxonomy_id,
                       self.__microbes_online_db,
                       self.__network_factories,
                       search_distances, scan_distances)


class OrganismBase:
    """The organism base class contains functionality that is likely to
    be the same among Organism implementations"""

    def __init__(self, code, network_factories):
        """Initialize the base class instance"""
        self.code = code
        self.__network_factories = network_factories
        self.__networks = None

    def networks(self):
        """returns this organism's networks"""
        if self.__networks == None:
            self.__networks = []
            for make_network in self.__network_factories:
                self.__networks.append(make_network(self))
        return self.__networks

    def thesaurus(self):
        """Returns a map containing the alias -> gene mappings"""
        raise Exception("please implement me")

    def feature_ids_for(self, gene_aliases):
        """Helper method to retrieve a list of feature_ids for the
        specified alias list"""
        synonyms = self.thesaurus()
        return [synonyms[alias] for alias in gene_aliases if alias in synonyms]


class Microbe(OrganismBase):
    """Abstraction of a microbial organism in cMonkey. It captures all
    organism-specific aspects. For now, we assume microbes only, but
    keep the interface generic so the algorithm will work on any type
    of organism"""

    # pylint: disable-msg=R0913,R0902
    def __init__(self, code, kegg_organism, rsat_info,
                 go_taxonomy_id, microbes_online_db,
                 network_factories,
                 search_distances, scan_distances):
        """create an Organism instance"""
        OrganismBase.__init__(self, code, network_factories)
        self.kegg_organism = kegg_organism
        self.__rsat_info = rsat_info
        self.__microbes_online_db = microbes_online_db
        self.go_taxonomy_id = go_taxonomy_id
        self.__search_distances = search_distances
        self.__scan_distances = scan_distances

        self.__synonyms = None  # lazy loaded
        self.__operon_mappings = None  # lazy loaded

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

    def features_for_genes(self, gene_aliases):
        """returns a map of features for the specified list of genes aliases
        used for operon information"""
        return util.ThesaurusBasedMap(
            self.thesaurus(),
            self.__read_features(self.feature_ids_for(gene_aliases)))

    def sequences_for_genes_search(self, gene_aliases, seqtype='upstream'):
        """The default sequence retrieval for microbes is to
        fetch their operon sequences"""
        return self.__operon_shifted_seqs_for(gene_aliases,
                                              self.__search_distances[seqtype])

    def sequences_for_genes_scan(self, gene_aliases, seqtype='upstream'):
        """The default sequence retrieval for microbes is to
        fetch their operon sequences"""
        return self.__operon_shifted_seqs_for(gene_aliases,
                                              self.__scan_distances[seqtype])

    def __operon_shifted_seqs_for(self, gene_aliases, distance):
        """returns a map of the gene_aliases to the feature-
        sequence tuple that they are actually mapped to.
        """
        def do_operon_shift():
            """Extract the (gene, head) pairs that are actually used"""
            operon_map = self.__operon_map()
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
            features = self.__read_features(unique_feature_ids)
            return self.__read_sequences(features, distance,
                                         st.extract_upstream)

        shifted_pairs = do_operon_shift()
        unique_seqs = unique_sequences(shifted_pairs)
        outseqs = {}
        for gene, head in shifted_pairs:
            outseqs[gene] = unique_seqs[head]
        return outseqs

    def __operon_map(self):
        """Returns the operon map for this particular organism.
        Microbes Online works on VNG names, but RSAT is working on
        feature ids, so this function also maps VNG names to feature ids"""
        if not self.__operon_mappings:
            pairs = mo.get_operon_pairs(self.__microbes_online_db, self)
            synonyms = self.thesaurus()
            self.__operon_mappings = {}
            for head, gene in pairs:
                self.__operon_mappings[synonyms[gene]] = synonyms[head]
        return self.__operon_mappings

    def thesaurus(self):
        """reads the thesaurus from a feature_names file. The thesaurus
        is also cached, because it is used many times
        """
        if not self.__synonyms:
            feature_names_dfile = util.DelimitedFile.create_from_text(
                self.__rsatdb().get_feature_names(self.species()),
                comment='--')
            self.__synonyms = thesaurus.create_from_rsat_feature_names(
                feature_names_dfile, [thesaurus.strip_vng_modification])
        return self.__synonyms

    def __read_features(self, feature_ids):
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
        dfile = util.DelimitedFile.create_from_text(
            self.__rsatdb().get_features(self.species()), comment='--')
        for line in dfile.lines():
            feature_id = line[0]
            if feature_id in feature_ids:
                features[feature_id] = read_feature(line)
        return features

    def __rsatdb(self):
        """internal method to return the RSAT db link"""
        return self.__rsat_info.rsatdb

    def __read_sequences(self, features, distance, extractor):
        """for each feature, extract and set its sequence"""
        def unique_contigs():
            """extract the unique contigs from the input features"""
            result = []
            for feature in features.values():
                if feature.location().contig not in result:
                    result.append(feature.location().contig)
            return result

        contig_seqs = {}
        sequences = {}
        for contig in unique_contigs():
            contig_seqs[contig] = self.__rsatdb().get_contig_sequence(
                self.species(), contig)

        for feature in features.values():
            location = feature.location()
            sequences[feature.id()] = extractor(
                contig_seqs[location.contig], location, distance)
        return sequences

    def __str__(self):
        result = "Organism Type: %s\n" % self.__class__.__name__
        result += (("Code: '%s'\nKEGG: '%s'\nRSAT: '%s'\nCOG: '%s'\n" +
                   "GO Taxonomy Id: %s\n") %
                   (self.code, self.kegg_organism, self.__rsat_info.species,
                    self.cog_organism(), self.go_taxonomy_id))
        return result


class GenericOrganism(OrganismBase):
    """Implementation of a generic organism that retrieves all of its data
    through data files"""

    def __init__(self, species_code,
                 thesaurus_filename, nw_factories,
                 seq_filenames,
                 search_distances, scan_distances):
        """Creates the organism"""
        OrganismBase.__init__(self, species_code, nw_factories)
        self.__seq_filenames = seq_filenames
        self.__thesaurus_filename = thesaurus_filename
        self.__search_distances = search_distances
        self.__scan_distances = scan_distances

        # lazy-loaded values
        self.__synonyms = None
        self.__seqs = {}

    def species(self):
        """Retrieves the species of this object"""
        return self.code

    def is_eukaryote(self):
        """Determines whether this object is an eukaryote"""
        return False

    def sequences_for_genes_search(self, gene_aliases, seqtype):
        """retrieve the sequences for the specified"""
        distance = self.__search_distances[seqtype]
        return self.__sequences_for_genes(seqtype, gene_aliases, distance)

    def sequences_for_genes_scan(self, gene_aliases, seqtype):
        """retrieve the sequences for the specified"""
        distance = self.__scan_distances[seqtype]
        return self.__sequences_for_genes(seqtype, gene_aliases, distance)

    def __sequences_for_genes(self, seqtype, gene_aliases, distance):
        """retrieves the specified sequences from the supplied genomic data"""
        if not seqtype in self.__seqs:
            logging.info('loading %s sequences' %seqtype)
            dfile = util.DelimitedFile.read(self.__seq_filenames[seqtype], sep=',')
            self.__seqs[seqtype] = {}
            for line in dfile.lines():
                self.__seqs[seqtype][line[0].upper()] = line[1].upper()
            logging.info('loaded %i %s sequences' \
                          %( len(self.__seqs[seqtype]), seqtype))
        result = {}
        for alias in gene_aliases:
            if alias in self.thesaurus():
                gene = self.thesaurus()[alias]
                if gene in self.__seqs[seqtype]:
                    result[gene] = self.__seqs[seqtype][gene]
                else:
                    #logging.warn("Gene '%s' not found in 3' UTRs", gene)
                    pass
            else:
                #logging.warn("Alias '%s' not in thesaurus !", alias)
                pass
        return result

    def thesaurus(self):
        """Reads the synonyms from the provided CSV file"""
        if not self.__synonyms:
            self.__synonyms = thesaurus.create_from_delimited_file2(
                self.__thesaurus_filename)
        return self.__synonyms


__all__ = ['make_kegg_code_mapper', 'make_go_taxonomy_mapper',
           'make_rsat_organism_mapper',
           'Organism', 'OrganismFactory']
