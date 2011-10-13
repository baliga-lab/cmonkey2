"""organism.py - organism-specific functionality in cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import logging
import util


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
        logging.info('RSAT - get_directory()')
        cache_file = "/".join([self.cache_dir, 'rsat_dir.html'])
        return util.read_url_cached("/".join([self.base_url,
                                              RsatDatabase.DIR_PATH]),
                                    cache_file)

    def get_organism(self, organism):
        """returns the file contents for the specified organism"""
        logging.info('RSAT - get_organism(%s)', organism)
        cache_file = "/".join([self.cache_dir, 'rsat_' + organism])
        return util.read_url_cached(
            "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                      RsatDatabase.ORGANISM_PATH]), cache_file)

    def get_organism_names(self, organism):
        """returns the specified organism name file contents"""
        logging.info('RSAT - get_organism_names(%s)', organism)
        cache_file = "/".join([self.cache_dir, 'rsatnames_' + organism])
        return util.read_url_cached(
            "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                      RsatDatabase.ORGANISM_NAMES_PATH]), cache_file)

    def get_ensembl_organism_names(self, organism):
        """returns the specified organism name file contents, using
        the EnsEMBL path"""
        logging.info('RSAT - get_ensembl_organism_names(%s)', organism)
        return util.read_url("/".join([self.base_url, RsatDatabase.DIR_PATH,
                                       organism + '_EnsEMBL',
                                       RsatDatabase.ORGANISM_NAMES_PATH]))

    def get_features(self, organism):
        """returns the specified organism's feature file contents
        Note: the current version only tries to read from feature.tab
        while the original cMonkey will fall back to cds.tab
        if that fails
        """
        #logging.info('RSAT - get_features(%s)', organism)
        cache_file = "/".join([self.cache_dir, organism + '_features'])
        return util.read_url_cached("/".join([self.base_url,
                                              RsatDatabase.DIR_PATH,
                                              organism,
                                              RsatDatabase.FEATURE_PATH]),
                               cache_file)

    def get_feature_names(self, organism):
        """returns the specified organism's feature name file contents"""
        #logging.info('RSAT - get_feature_names(%s)', organism)
        cache_file = "/".join([self.cache_dir, organism + '_feature_names'])
        return util.read_url_cached(
            "/".join([self.base_url,
                      RsatDatabase.DIR_PATH,
                      organism,
                      RsatDatabase.FEATURE_NAMES_PATH]),
            cache_file)

    def get_contig_sequence(self, organism, contig):
        """returns the specified contig sequence"""
        #logging.info('RSAT - get_contig_sequence(%s, %s)',
        #             organism, contig)
        cache_file = "/".join([self.cache_dir, organism + '_' + contig])
        url = "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                        'genome', contig + '.raw'])
        return util.read_url_cached(url, cache_file).upper()


__all__ = ['RsatDatabase']
