# vi: sw=4 ts=4 et:
"""organism.py - organism-specific functionality in cMonkey
This module captures a microbial organism that receives data
from Microbes Online and RSAT

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import logging
import collections

import cmonkey.thesaurus as thesaurus
import cmonkey.util as util
import cmonkey.seqtools as st
import cmonkey.microbes_online as mo
import cmonkey.patches as patches

# requires biopython
from Bio import SeqIO


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
                 ratios=None, synonyms=None, fasta_file=None):
        """create an Organism instance"""
        # microbe-specific network factories need access to synonyms
        # and rsat info, so initialize them here before the base class
        # init
        self.__synonyms = synonyms
        self.__rsat_info = rsat_info
        OrganismBase.__init__(self, code, network_factories, ratios=ratios)
        self.kegg_organism = kegg_organism
        self.go_taxonomy_id = go_taxonomy_id
        self.search_distances = search_distances
        self.scan_distances = scan_distances
        self.fasta_file = fasta_file
        if self.fasta_file is not None:
            self.sequence_source = FASTASequenceSource(self, self.fasta_file)
        else:
            self.sequence_source = RSATOrganismSequenceSource(self)

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
                                          int(line[4].lstrip('<>')),
                                          int(line[5].lstrip('<>')),
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

    def sequences_for_genes_search(self, genes, seqtype='upstream'):
        """The default sequence retrieval for microbes is to
        fetch their operon sequences"""
        return self.sequence_source.seqs_for(genes, self.search_distances[seqtype])

    def sequences_for_genes_scan(self, genes, seqtype='upstream'):
        """The default sequence retrieval for microbes is to
        fetch their operon sequences"""
        return self.sequence_source.seqs_for(genes, self.scan_distances[seqtype])


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
                 use_operons=True, ratios=None, synonyms=None, fasta_file=None):
        """create an Organism instance"""
        RSATOrganism.__init__(self, code, kegg_organism,
                              rsat_info, go_taxonomy_id, network_factories,
                              search_distances, scan_distances, ratios, synonyms,
                              fasta_file)
        self.use_operons = use_operons
        self.__microbes_online_db = microbes_online_db
        self.__operon_mappings = None  # lazy loaded

    def operon_map(self):
        """Returns the operon map for this particular organism.
        Microbes Online works on VNG names, but RSAT is working on
        feature ids, so this function also maps VNG names to feature ids"""
        if not self.__operon_mappings:
            pairs = mo.get_operon_pairs(self.__microbes_online_db, self)
            synonyms = self.thesaurus()
            self.__operon_mappings = {synonyms[gene]: synonyms[head] for head, gene in pairs}
        return self.__operon_mappings


class FASTASequenceSource:
    """FASTA file based sequence source"""

    def __init__(self, organism, filepath):
        self.organism = organism
        self.seqmap = None
        with open(filepath) as infile:
            self.fasta_records = [r for r in SeqIO.parse(infile, 'fasta')]

    def seqs_for(self, gene_aliases, distances):
        def seq2str(seq):
            return str(seq.upper()).replace('X', 'N')

        synonyms = self.organism.thesaurus()
        if self.seqmap is None:
            self.seqmap = {synonyms[r.id]: r for r in self.fasta_records if r.id in synonyms}
        result = {}
        for gene in gene_aliases:
            if gene in synonyms and synonyms[gene] in self.seqmap:
                fasta_record = self.seqmap[synonyms[gene]]
                result[synonyms[gene]] = (st.Location('', 0, 0, False),
                                          seq2str(fasta_record.seq.upper()))
            else:
                logging.warn("'%s' not in the sequences !!!!", gene)
        return result


class RSATOrganismSequenceSource:
    """Default sequence source for Microbes"""

    def __init__(self, organism):
        self.organism = organism

    def seqs_for(self, gene_aliases, distance):
        """returns a map of the gene_aliases to the feature-
        sequence tuple that they are actually mapped to.
        """
        def do_operon_shift():
            """Extract the (gene, head) pairs that are actually used"""
            operon_map = self.organism.operon_map()
            synonyms = self.organism.thesaurus()
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
            features = self.organism.read_features(unique_feature_ids)
            return self.organism.read_sequences(features, distance,
                                                st.extract_upstream)

        if self.organism.use_operons:
            shifted_pairs = do_operon_shift()
        else:
            # if operons should not be used, we simply map
            # the gene heads to themselves
            synonyms = self.organism.thesaurus()
            valid_genes = [synonyms[alias]
                           for alias in gene_aliases if alias in synonyms]
            shifted_pairs = [(gene, gene) for gene in valid_genes]

        unique_seqs = unique_sequences(shifted_pairs)
        result = {}
        skipped = set()
        for gene, head in shifted_pairs:
            if head in unique_seqs:
                result[gene] = unique_seqs[head]
            else:
                skipped.add(head)
        #return {gene: unique_seqs[head] for gene, head in shifted_pairs}
        if len(skipped) > 0:
            logging.warn("%d genes were skipped while reading their sequences from the RSAT files", len(skipped))
        return result


__all__ = ['RSATOrganism', 'Microbe']
