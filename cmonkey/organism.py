"""organism.py - organism-specific functionality in cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
from util import DelimitedFileMapper


class KeggCodeMapper:
    """A class to map an organism code to a KEGG organism name. It uses
    a KEGG taxonomy file to do this"""

    def __init__(self, delimited_file):
        """Creates an instance of the mapper class using a DelimitedFile"""
        self.file_mapper = DelimitedFileMapper(delimited_file, 1, 3)

    def get_organism(self, code):
        """lookup the organism code to return the full name"""
        return self.file_mapper.lookup(code)


class GoTaxonomyMapper:
    """A class to map an organism name to a GO taxonomy id"""

    def __init__(self, delimited_file):
        """Creates an instance of the mapper class"""
        self.file_mapper = DelimitedFileMapper(delimited_file, 0, 1)

    def get_taxonomy_id(self, organism_name):
        """look up the taxonomy id for a given organism name. Note that
        organism names are not necessarily unique.
        This will return the first one found"""
        return self.file_mapper.lookup(organism_name)


class OrganismFactory:
    """Factory to create an organism. Construction of an organism
    instance is relatively complex and costly, so it is coordinated
    here. Information has to be pulled together from various databases
    which are provided to the factory as configuration parameters.
    """

    def __init__(self, kegg_organism_mapper,
                 rsat_organism_mapper,
                 go_taxonomy_mapper):
        self.kegg_organism_mapper = kegg_organism_mapper
        self.rsat_organism_mapper = rsat_organism_mapper
        self.go_taxonomy_mapper = go_taxonomy_mapper

    def create(self, organism_code):
        """factory method to create an organism from a code"""
        kegg_organism = self.kegg_organism_mapper.get_organism(organism_code)
        rsat_organism = self.rsat_organism_mapper.get_organism(kegg_organism)
        go_taxonomy_id = self.go_taxonomy_mapper.get_taxonomy_id(rsat_organism)
        return Organism(organism_code, kegg_organism, rsat_organism,
                        go_taxonomy_id)


class Organism:
    """Abstraction of an organism in cMonkey. It captures all organism-specific
    aspects"""
    def __init__(self, code, kegg_organism, rsat_organism, go_taxonomy_id):
        """create an Organism instance"""
        self.code = code
        self.kegg_organism = kegg_organism
        self.rsat_organism = rsat_organism
        self.go_taxonomy_id = go_taxonomy_id

    def is_prokaryote(self):
        """determine whether this organism is a prokaryote"""
        return True

    def is_eukaryote(self):
        """determine whether this organism is a eukaryote"""
        return not self.is_prokaryote()


__all__ = ['KeggCodeMaper', 'GoTaxonomyMapper', 'Organism', 'OrganismFactory']
