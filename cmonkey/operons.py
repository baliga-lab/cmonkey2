"""operons.py - operon prediction loading cMonkey
For organisms that have operons, cMonkey can read use operon predictions
to improve its scoring. This module defines access to prediction data through
various services.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
from util import read_url_cached, DelimitedFile


MICROBES_ONLINE_BASE_URL = 'http://www.microbesonline.org'


class MicrobesOnline:
    """Interface to Microbes Online web service"""

    def __init__(self, base_url=MICROBES_ONLINE_BASE_URL,
                 cache_dir='cache'):
        """creates a MicrobesOnline service instance"""
        self.base_url = base_url
        self.cache_dir = cache_dir

    def get_operon_predictions_for(self, organism_id):
        """Retrieve operon predictions for the specified organism"""
        url = '/'.join([self.base_url, 'operons',
                       'gnc%s.named' % str(organism_id)])
        cache_file = '/'.join([self.cache_dir,
                              'gnc%s.named' % str(organism_id)])
        return read_url_cached(url, cache_file)


def read_from_microbes_online(microbes_online, organism_id):
    """reads operon predictions from Microbes Online"""
    preds = microbes_online.get_operon_predictions_for(organism_id)
    dfile = DelimitedFile.create_from_text(preds, has_header=True)
    return [(line[2], line[3]) for line in dfile.lines() if line[6] == 'TRUE']


__all__ = ['MicrobesOnline', 'read_from_microbes_online']
