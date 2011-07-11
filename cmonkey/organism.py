"""organism.py - organism-specific functionality in cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
from util import DelimitedFileMapper, best_matching_links
import urllib
import re


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
    """A class to map an RSAT organism name to a GO taxonomy id"""

    def __init__(self, delimited_file):
        """Creates an instance of the mapper class"""
        self.file_mapper = DelimitedFileMapper(delimited_file, 0, 1)

    def get_taxonomy_id(self, organism_name):
        """look up the taxonomy id for a given organism name. Note that
        organism names are not necessarily unique.
        This will return the first one found"""
        return self.file_mapper.lookup(organism_name)


class RsatOrganismMapper:
    """A class to map a KEGG organism name to an RSAT organism name"""

    def __init__(self, rsat_database):
        """Creates a mapper instance with an attached RSAT database"""
        self.rsatdb = rsat_database

    def get_organism(self, kegg_organism):
        """find the RSAT organism which best matches the given KEGG
        organism name"""
        return best_matching_links(
            kegg_organism,
            self.rsatdb.get_directory_html())[0].rstrip('/')

    def is_eukaryotic(self, rsat_organism):
        """determine whether the given RSAT organism is eukaryotic"""
        organism_text = self.rsatdb.get_organism_file(rsat_organism)
        return re.search('Eukaryota', organism_text)


class OrganismFactory:
    """Factory to create an organism. Construction of an organism
    instance is relatively complex and costly, so it is coordinated
    here. Information has to be pulled together from various databases
    which are provided to the factory as configuration parameters.
    """

    def __init__(self, kegg_organism_mapper,
                 rsat_organism_mapper,
                 go_taxonomy_mapper):
        """create a OrganismFactory instance"""
        self.kegg_organism_mapper = kegg_organism_mapper
        self.rsat_organism_mapper = rsat_organism_mapper
        self.go_taxonomy_mapper = go_taxonomy_mapper

    def create(self, organism_code):
        """factory method to create an organism from a code"""
        kegg_organism = self.kegg_organism_mapper.get_organism(organism_code)
        rsat_organism = self.rsat_organism_mapper.get_organism(kegg_organism)
        go_taxonomy_id = self.go_taxonomy_mapper.get_taxonomy_id(rsat_organism)
        is_eukaryote = self.rsat_organism_mapper.is_eukaryote(rsat_organism)
        if is_eukaryote:
            return Eukaryote(organism_code, kegg_organism, rsat_organism,
                             go_taxonomy_id)
        else:
            return Prokaryote(organism_code, kegg_organism, rsat_organism,
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

    def get_cog_organism(self):
        """returns the COG organism name"""
        return self.code.capitalize()


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


class RsatDatabase:
    """abstract interface to access an RSAT mirror"""
    DIR_PATH = 'data/genomes'
    ORGANISM_FILE_PATH = 'genome/organism.tab'

    def __init__(self, base_url):
        """create an RsatDatabase instance based on a mirror URL"""
        self.base_url = base_url

    def get_directory_html(self):
        """returns the HTML page for the directory listing"""
        return urllib.urlopen("/".join([self.base_url,
                                        RsatDatabase.DIR_PATH])).read()

    def get_organism_file(self, organism):
        """returns the file contents for the specified organism"""
        return urllib.urlopen(
            "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                      RsatDatabase.ORGANISM_FILE_PATH])).read()


__all__ = ['KeggCodeMaper', 'GoTaxonomyMapper', 'Organism', 'OrganismFactory',
           'RsatDatabase']
