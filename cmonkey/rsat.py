"""organism.py - organism-specific functionality in cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import urllib


class DocumentNotFound(Exception):
    """An exception indicating that the requested document does not exist"""
    pass


class CMonkeyURLopener(urllib.FancyURLopener):
    """An URL opener that can detect 404 errors"""

    def http_error_default(self, url, fp, errcode, errmsg, headers):
        # pylint: disable-msg=R0913
        # pylint: disable-msg=C0103
        """overriding the default error handling method to handle HTTP 404
        errors"""
        if (errcode == 404):
            raise DocumentNotFound(url)

        # call super class handler.
        # note that urllib.FancyURLopener is not a new-style class
        return urllib.FancyURLopener.http_error_default(
            self, url, fp, errcode, errmsg, headers)


def read_url(url):
    """convenience method to read a document from a URL using the
    CMonkeyURLopener"""
    return CMonkeyURLopener().open(url).read()


class RsatDatabase:
    """abstract interface to access an RSAT mirror"""
    DIR_PATH = 'data/genomes'
    ORGANISM_PATH = 'genome/organism.tab'
    ORGANISM_NAMES_PATH = 'genome/organism_names.tab'
    FEATURE_PATH = 'genome/feature.tab'

    def __init__(self, base_url):
        """create an RsatDatabase instance based on a mirror URL"""
        self.base_url = base_url

    def get_directory(self):
        """returns the HTML page for the directory listing"""
        return read_url("/".join([self.base_url,
                                        RsatDatabase.DIR_PATH]))

    def get_organism(self, organism):
        """returns the file contents for the specified organism"""
        return read_url(
            "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                      RsatDatabase.ORGANISM_PATH]))

    def get_organism_names(self, organism):
        """returns the specified organism name file contents"""
        return read_url(
            "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                      RsatDatabase.ORGANISM_NAMES_PATH]))

    def get_ensembl_organism_names(self, organism):
        """returns the specified organism name file contents, using
        the EnsEMBL path"""
        return read_url("/".join([self.base_url, RsatDatabase.DIR_PATH,
                                  organism + '_EnsEMBL',
                                  RsatDatabase.ORGANISM_NAMES_PATH]))

    def get_features(self, organism):
        """returns the specified organism's feature file contents"""
        return read_url("/".join([self.base_url, RsatDatabase.DIR_PATH,
                                  organism,
                                  RsatDatabase.FEATURE_PATH]))

__all__ = ['RsatDatabase']
