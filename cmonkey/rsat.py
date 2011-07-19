"""organism.py - organism-specific functionality in cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import os
from util import read_url, read_url_cached, CMonkeyURLopener


class RsatDatabase:
    """abstract interface to access an RSAT mirror"""
    DIR_PATH = 'data/genomes'
    ORGANISM_PATH = 'genome/organism.tab'
    ORGANISM_NAMES_PATH = 'genome/organism_names.tab'
    FEATURE_PATH = 'genome/feature.tab'
    FEATURE_NAMES_PATH = 'genome/feature_names.tab'

    def __init__(self, base_url, cache_dir):
        """create an RsatDatabase instance based on a mirror URL"""
        self.base_url = base_url
        self.cache_dir = cache_dir.rstrip('/')

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
        """returns the specified organism's feature file contents
        Note: the current version only tries to read from feature.tab
        while the original cMonkey will fall back to cds.tab
        if that fails
        """
        return read_url("/".join([self.base_url, RsatDatabase.DIR_PATH,
                                  organism,
                                  RsatDatabase.FEATURE_PATH]))

    def get_feature_names(self, organism):
        """returns the specified organism's feature name file contents"""
        return read_url("/".join([self.base_url, RsatDatabase.DIR_PATH,
                                  organism,
                                  RsatDatabase.FEATURE_NAMES_PATH]))

    def get_contig_sequence(self, organism, contig):
        """returns the specified contig sequence"""
        cache_file = "/".join([self.cache_dir, organism + '_' + contig])
        url = "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                        'genome', contig + '.raw'])
        return read_url_cached(url, cache_file).upper()

    def cache_contig_sequence(self, organism, contig):
        """downloads the specified contig sequence to the cache directory
        if it does not yet exist"""
        cache_file = "/".join([self.cache_dir, organism + '_' + contig])
        if not os.path.exists(cache_file):
            url = "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                            'genome', contig + '.raw'])
            print "RETRIEVE URL: %s" % url
            CMonkeyURLopener().retrieve(url, cache_file)
        else:
            print "found existing cache file: '%s'" % cache_file

__all__ = ['RsatDatabase']
