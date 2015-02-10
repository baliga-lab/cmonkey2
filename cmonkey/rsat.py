# vi: sw=4 ts=4 et:
"""organism.py - organism-specific functionality in cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import logging
import util
import StringIO
import re
import patches
import os

class RsatFiles:
    """This class implements the same service functions as RsatDatabase, but
    takes the data from files"""
    def __init__(self, dirname, basename, taxonomy_id, feature_name, url):
        self.dirname = dirname
        self.taxonomy_id = taxonomy_id
        self.basename = basename
        self.feature_name = feature_name
        self.url = url

    def get_taxonomy_id(self, organism):
        return self.taxonomy_id

    def get_rsat_organism(self, kegg_organism):
        return self.basename

    def get_rsat_getURL(self):
        return self.url

    def get_rsat_featureName(self):
        return self.featureName

    def get_features(self, organism, original=True):
        if original:
            #path = os.path.join(self.dirname, 'features.tab')
            path = os.path.join(self.dirname, self.feature_name + '.tab')
        else:
            path = os.path.join(self.dirname, organism + '_' + self.feature_name)
        with open(path) as infile:
            return infile.read()

    def get_feature_names(self, organism, original=True):
        import pdb
        pdb.set_trace()
        if original:
            path = os.path.join(self.dirname, self.feature_name + '_names.tab')
        else:
            path = os.path.join(self.dirname, organism + '_' + self.feature_name + '_names')
        with open(path) as infile:
            return infile.read()

    def get_contig_sequence(self, organism, contig, original=True):
        if original:
            path = os.path.join(self.dirname, contig + '.tab')
        else:
            path = os.path.join(self.dirname, organism + '_' + contig)
        with open(path) as infile:
            seqstr = infile.read().upper()
            return join_contig_sequence(seqstr)


class RsatDatabase:
    """abstract interface to access an RSAT mirror"""
    DIR_PATH = 'data/genomes'
    ORGANISM_PATH = 'genome/organism.tab'
    ORGANISM_NAMES_PATH = 'genome/organism_names.tab'
    #FEATURE_PATH = 'genome/feature.tab'
    #FEATURE_NAMES_PATH = 'genome/feature_names.tab'

    def __init__(self, base_url, cache_dir, kegg_species, ncbi_code, feature_name='feature'):
    
        """create an RsatDatabase instance based on a mirror URL"""
        self.base_url = base_url
        self.cache_dir = cache_dir.rstrip('/')
        self.kegg_species = kegg_species
        self.ncbi_code = ncbi_code
        self.feature_name = feature_name
        self.feature_path = 'genome/' + feature_name + '.tab'
        self.feature_names_path = 'genome/' + feature_name + '_names.tab'

    def __get_ncbi_code(self, rsat_organism):
        """retrieve NCBI code from organism.tab file"""
        cache_file = "/".join([self.cache_dir, '%s.tab' % rsat_organism])
        text = util.read_url_cached("/".join([self.base_url,
                                              RsatDatabase.DIR_PATH,
                                              rsat_organism,
                                              RsatDatabase.ORGANISM_PATH]),
                                    cache_file)
        spec = [line for line in text.split('\n') if not line.startswith('--')][0]
        return spec.strip().split('\t')[0]

    def get_rsat_organism(self, kegg_organism):
        """returns the HTML page for the directory listing"""
        logging.info('RSAT - get_directory()')
        cache_file = "/".join([self.cache_dir, 'rsat_dir.html'])
        text = util.read_url_cached("/".join([self.base_url,
                                              RsatDatabase.DIR_PATH]),
                                    cache_file)
        suggestion1 = util.best_matching_links(self.kegg_species, text)[0].rstrip('/')
        suggestion2 = util.best_matching_links(kegg_organism, text)[0].rstrip('/')
        if suggestion1 != suggestion2:
            ncbi_code1 = self.__get_ncbi_code(suggestion1)
            ncbi_code2 = self.__get_ncbi_code(suggestion1)
            if ncbi_code1 == self.ncbi_code:
                return suggestion1
            elif ncbi_code2 == self.ncbi_code:
                return suggestion2
            else:
                logging.warn("can't find the correct RSAT mapping !")
                return suggestion1
        else:
            ncbi_code = self.__get_ncbi_code(suggestion1)
            if ncbi_code == self.ncbi_code:
                return suggestion1
            else:
                logging.warn("can't find the correct RSAT mapping !")
                return suggestion1

    def get_taxonomy_id(self, organism):
        """returns the specified organism name file contents"""
        logging.info('RSAT - get_organism_names(%s)', organism)
        cache_file = "/".join([self.cache_dir, 'rsatnames_' + organism])
        text = util.read_url_cached(
            "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                      RsatDatabase.ORGANISM_NAMES_PATH]), cache_file)
        organism_names_dfile = util.dfile_from_text(text, comment='--')
        return patches.patch_ncbi_taxonomy(organism_names_dfile.lines[0][0])

    def get_features(self, organism):
        """returns the specified organism's feature file contents
        Note: the current version only tries to read from feature.tab
        while the original cMonkey will fall back to cds.tab
        if that fails
        """
        logging.info('RSAT - get_features(%s)', organism)
        cache_file = "/".join([self.cache_dir, organism + '_' + self.feature_name])
        uCache = util.read_url_cached("/".join([self.base_url, RsatDatabase.DIR_PATH, organism, self.feature_path]), cache_file)

        #Make sure that the fields are in the correct order
        #Later parts assume that the features file will have the following columns
        fieldOrder = ['id', 'type', 'name', 'contig', 'start_pos', 'end_pos', 'strand']

        uCache = uCache.split('\n')
        #Remove any blank lines
        while "" in uCache:
            uCache.remove("")

        idxs = {} #Dictionary to store field idxs
        targIdx = [] #The ordered list of columns for output
        outString = "" #This will be the new data
        for line in uCache:
            try:
                line = line + '\n'
            except:
                continue
        
            lineParts = line.split()
            if lineParts[0] == '--':
                if lineParts[1] == 'field':
                        idxs[lineParts[3]] = lineParts[2]
                        if lineParts[3] in fieldOrder:
                                newIdx = str(fieldOrder.index(lineParts[3]) + 1)
                                outString = outString + lineParts[0] + " " + lineParts[1] + " " + newIdx + '\t' + lineParts[3] + '\n'
                else:
                        outString = outString + line
            else:
                if (len(targIdx) == 0):
                        #Create the targIdx
                        for curField in fieldOrder:
                                targIdx.append(int(idxs[curField])-1)
                outline = ""
                for curTarg in targIdx:
                        outline = outline + lineParts[curTarg] + '\t'
                #Some RSAT files have a contig with ':'s instead of '_'s
                outline = outline.replace(':','_')
                #Now strip trailing \t
                outline = ''.join(outline.rsplit('\t', 1))
                outString = outString + outline + '\n'

        #To Do: Overwrite cache file & add early check to see if we need the sub
        return outString

    def get_feature_names(self, organism):
        """returns the specified organism's feature name file contents"""
        #logging.info('RSAT - get_feature_names(%s)', organism)
        cache_file = "/".join([self.cache_dir, organism + '_' + self.feature_name + '_names'])
        uCach = util.read_url_cached(
            "/".join([self.base_url,
                      RsatDatabase.DIR_PATH,
                      organism,
                      self.feature_names_path]),
            cache_file)
        return uCach

    def get_contig_sequence(self, organism, contig):
        """returns the specified contig sequence"""
        logging.info('RSAT - get_contig_sequence(%s, %s)',
                    organism, contig)
        cache_file = "/".join([self.cache_dir, organism + '_' + contig])
        url = "/".join([self.base_url, RsatDatabase.DIR_PATH, organism,
                        'genome', contig + '.raw'])
    
        #10-07-14 Crashed here with URL timeout.  Maybe RSAT limits downloads?
        #  On 10-08-14 I could download the other files with pdb.set_trace()
        #  Maybe all I will need is a pause between files?
        try:
            seqstr = util.read_url_cached(url, cache_file).upper()
        except:
            print "Error downloading file: " + url
            print "RSAT occasionally has connectivity problems."
            print "Try again later, or try a different RSAT mirror"
            print "useing the parameter --rsat_base_url"
        return join_contig_sequence(seqstr)

def join_contig_sequence(seqstr):
    """we take the safer route and assume that the input could
    be separated out into lines"""
    buf = StringIO.StringIO(seqstr)
    result = ''
    for line in buf:
        result += line.strip()
    return result

__all__ = ['RsatDatabase']
