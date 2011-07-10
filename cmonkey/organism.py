"""organism.py - organism-specific functionality in cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


def get_kegg_organism_for_code(kegg_taxonomy_file, code):
    """using a KEGG taxonomy file, lookup the organism code to
    return the full name. taxonomy_file should be an instance of
    DelimitedFile"""
    for line in kegg_taxonomy_file.get_lines():
        if line[1] == code:
            return line[3]
    return None


def get_go_taxonomy_id(go_taxonomy_file, organism_name):
    """using a GO proteome2taxid file, look up the taxonomy id for a given
    organism name. Note that organism names are not necessarily unique.
    This will return the first one found"""
    for line in go_taxonomy_file.get_lines():
        if line[0] == organism_name:
            return line[1]
    return None


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


__all__ = ['get_organism_for_code', 'get_go_taxonomy_id', 'Organism',
           'OrganismFactory']
