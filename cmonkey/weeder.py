# vi: sw=4 ts=4 et:
"""weeder.py - cMonkey weeder interface module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import subprocess as sp
import logging
import re
import os
import cmonkey.pssm as pssm
import cmonkey.weederlauncher as weederlauncher

# Python2/Python3 compatibility
try:
    xrange
except NameError:
    xrange = range


LAUNCHER = 'weederlauncher'
LLR_VALUE = 'NA'


class Site:
    """A site value object holding entries from the wee file"""
    def __init__(self, gene, reverse, site, start, match):
        """create Site object"""
        self.gene = gene
        self.is_reverse = reverse
        self.site = site
        self.start = start
        self.match = match

    def __str__(self):
        return ("gene = %s, rev = %d, site = %s, start = %d, match = %s" %
                (self.gene, self.is_reverse, self.site, self.start,
                 self.match))

    def __repr__(self):
        return self.__str__()

def run_weeder(fasta_file, params, config_params, bgmodel):
    if not os.path.exists(fasta_file):
        logging.warning("Weeder FASTA file %s not found! Skipping")
        return []
    meme_outfile = '%s/meme-out-%04d-%04d.txt' % (params.outdir, params.iteration,
                                                  params.cluster)

    """run the weeder command and interpret its result"""
    def write_f1_file(pssm_num, apssm, num_sites):
        """writes the pssm to the .f1 file"""
        with open('%s.%d.f1' % (fasta_file, pssm_num), 'w') as outfile:
            outfile.write('%d,%s,%f,%d\n' % (apssm.sequence_length(),
                                             LLR_VALUE,
                                             apssm.e_value,
                                             num_sites))
            for row in xrange(apssm.sequence_length()):
                outfile.write('%f,%f,%f,%f\n' % (apssm[row][0],
                                                 apssm[row][1],
                                                 apssm[row][2],
                                                 apssm[row][3]))

    def write_f2_file(pssm_num, apssm):
        """writes the pssm to the .f2 file"""
        with open('%s.%d.f2' % (fasta_file, pssm_num), 'w') as outfile:
            outfile.write(',gene,strand,start,p.value,site\n')
            for index in xrange(len(apssm.sites)):
                site = apssm.sites[index]
                strand = '+'
                if site.is_reverse:
                    strand = '-'

                outfile.write('%d,%s,%s,%d,%f,%s\n' %
                              (index + 1, site.gene, strand, site.start,
                               site.match, site.site))

    def write_meme_file(pssms):
        """writes a PSSM file to be read by meme"""
        with open(meme_outfile, 'w') as outfile:
            outfile.write("""MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A %.3f C %.3f G %.3f T %.3f

""" % (bgmodel[0]['A'], bgmodel[0]['C'], bgmodel[0]['G'], bgmodel[0]['T']))
            for motif_num, apssm in enumerate(pssms):
                outfile.write('MOTIF m%d\n' % motif_num)
                outfile.write(apssm.to_logodds_string())
                outfile.write('\n')

    __launch_weeder(fasta_file, config_params)
    pssms = [pssm8 for pssm8 in __read_pssms_for(fasta_file)
             if pssm8.sequence_length() == 8]
    for index in xrange(len(pssms)):
        unique_target_sites = []
        for site in pssms[index].sites:
            if site.gene not in unique_target_sites:
                unique_target_sites.append(site.gene)
        num_sites = len(unique_target_sites)
        write_f1_file(index + 1, pssms[index], num_sites)
        write_f2_file(index + 1, pssms[index])
    write_meme_file(pssms)
    return meme_outfile, pssms


def __launch_weeder(fasta_file, config_params):
    """launch weeder command"""
    if config_params['Weeder']['orgcode']:
        orgcode = config_params['Weeder']['orgcode']
    else:
        orgcode = config_params['organism_code'].upper()

    freqfile_dir = config_params['Weeder']['freqfile_dir']
    weederlauncher.run_small_analysis(fasta_file, orgcode, 50,
                                      reverse=True,
                                      multi=False,
                                      allseqs=False,
                                      ffdir=freqfile_dir)


def __read_pssms_for(fasta_file):
    """parses the output from weeder into PSSMS"""
    reader = WeederReader(fasta_file + ".wee", 'pssm')
    reader.read()
    return reader.pssms()


class WeederReader:
    """Reads PSSMs from .wee files"""
    def __init__(self, wee_file, base_name):
        """create a reader instance"""
        self.__wee_file = wee_file
        self.__lines = None
        self.__current_index = 0
        self.__top_hit6 = None
        self.__top_hit8 = None
        self.__sequence_names = None
        self.__pssms = None
        self.__base_name = base_name

    def top_hit6(self):
        """returns the top hit of length 6, a (sequence, score) pair"""
        return self.__top_hit6

    def top_hit8(self):
        """returns the top hit of length 8, a (sequence, score) pair"""
        return self.__top_hit8

    def sequence_names(self):
        """returns the sequence names read"""
        return self.__sequence_names

    def pssms(self):
        """returns the PSSMs read"""
        return self.__pssms

    def read(self):
        """takes a wee file name and parses the contents into a list"""
        with open(self.__wee_file) as infile:
            self.__lines = infile.readlines()
        self.__top_hit6 = self.__locate_top_entry_for_length(6)
        self.__top_hit8 = self.__locate_top_entry_for_length(8)
        self.__sequence_names = self.__read_sequence_names()
        self.__pssms = self.__read_pssms()

    def __current_line(self):
        """returns the line currently pointed to"""
        return self.__lines[self.__current_index]

    def __input_end_reached(self):
        """determines whether the end of input was reached"""
        return self.__current_index >= len(self.__lines)

    def __find_line_starting_with(self, start):
        """finds the next line that starts with a certain string"""
        while (not self.__input_end_reached() and
               not self.__current_line().startswith(start)):
            self.__current_index += 1
        if self.__input_end_reached():
            raise Exception("Line starting with '", start, "' not found")

    def __find_line_matching(self, regex):
        """finds the next line that starts with a certain string"""
        while (not self.__input_end_reached() and
               not re.match(regex, self.__current_line())):
            self.__current_index += 1
        if self.__input_end_reached():
            raise Exception("Line matching '", regex, "' not found")

    def __next_nonempty_line(self):
        """finds the next line that starts with a certain string"""
        self.__current_index += 1
        while (not self.__input_end_reached() and
               len(self.__current_line().strip()) == 0):
            self.__current_index += 1
        if self.__input_end_reached():
            raise Exception("no more non-empty line found")

    def __locate_top_entry_for_length(self, seqlen):
        """locate the top entry for the specified sequence length"""
        self.__find_line_starting_with('Searching for motifs of length ' +
                                       str(seqlen))
        return self.__top_hit_entry()

    def __top_hit_entry(self):
        """returns the top hit entry directly following the specified index"""
        self.__find_line_starting_with('1)')
        comps = self.__current_line().split(' ')
        return (comps[1], float(comps[2]))

    def __read_sequence_names(self):
        """reads and returns the sequence names"""
        result = []
        self.__find_line_starting_with('Your sequences:')
        while (not self.__input_end_reached() and
               not self.__current_line().startswith('**** MY ADVICE ****')):
            comps = self.__current_line().strip().split(' ')
            if len(comps) == 4 and not self.__current_line().startswith('****'):
                result.append(comps[3][1:])
            self.__current_index += 1
        return result

    def __read_pssms(self):
        """reads and returns the PSSM's"""
        result = []
        self.__find_line_starting_with(
            ' *** Interesting motifs (highest-ranking)')
        self.__next_nonempty_line()
        while (not self.__input_end_reached() and
               not re.match('.*not highest-ranking.*$',
                            self.__current_line())):
            pssm_name = self.__base_name + '_' + self.__current_line().strip()
            sites = self.__read_sites()
            matrix = self.__read_frequency_matrix()
            if len(matrix) == 6:
                top_hit = self.__top_hit6
            else:
                top_hit = self.__top_hit8
            result.append(pssm.Pssm(pssm_name, matrix, e_value=top_hit[1],
                                    sites=sites))
            self.__current_index += 2  # skip over the empty line and the ====
            self.__next_nonempty_line()  # and position to the next PSSM
        return result

    def __read_sites(self):
        """reads the sites from the best occurrences section"""
        result = []
        self.__find_line_starting_with('Best occurrences:')
        self.__current_index += 2
        while len(self.__current_line().strip()) > 0:
            comps = re.split('[ \t]+', self.__current_line().strip())
            result.append(Site(self.__sequence_names[int(comps[0]) - 1],
                               comps[1] != '+', comps[2], int(comps[3]),
                               float(comps[4].lstrip('(').rstrip(')'))))
            self.__current_index += 1
        return result

    def __read_frequency_matrix(self):
        """reads the frequency matrix for one pssm"""
        self.__find_line_matching('\s+Frequency Matrix.*$')
        self.__current_index += 3  # skips header infos
        matrix = []
        while len(self.__current_line().strip()) > 0:
            comps = re.split('[ \t]+', self.__current_line().strip())
            all_occurrences = [float(value) for value in comps[1:5]]
            col_sum = sum(all_occurrences)
            matrix.append([(occ / col_sum) for occ in all_occurrences])
            self.__current_index += 1
        return matrix
