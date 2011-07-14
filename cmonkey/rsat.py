"""organism.py - organism-specific functionality in cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import urllib


class RsatDatabase:
    """abstract interface to access an RSAT mirror"""
    DIR_PATH = 'data/genomes'
    ORGANISM_FILE_PATH = 'genome/organism.tab'
    ORGANISM_NAMES_FILE_PATH = 'genome/organism_names.tab'

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

    def get_organism_names_file(self, organism):
        """returns the specified organism name file contents"""
        return urllib.urlopen(
            "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                      RsatDatabase.ORGANISM_FILE_PATH])).read()

__all__ = ['RsatDatabase']
