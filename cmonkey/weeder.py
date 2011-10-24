"""weeder.py - cMonkey weeder interface module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import subprocess as sp
import logging
import re
import pssm

LAUNCHER = 'weederlauncher'


def run_weeder(fasta_file):
    """run the weeder command and interpret its result"""
    __launch_weeder(fasta_file)
    return __read_pssms_for(fasta_file)


def __launch_weeder(fasta_file):
    """launch weeder command"""
    command = [LAUNCHER, fasta_file, 'HS3P', 'small', 'T50']
    retcode = 1
    with open('weeder.log', 'w') as logfile:
        logging.info("running weeder on '%s'", fasta_file)
        weederproc = sp.Popen(command, stdout=logfile, stderr=sp.STDOUT)
        retcode = weederproc.wait()
        logging.info("Weeder finished, return code: %d", retcode)
    return retcode


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
        print "lines read: ", len(self.__lines)
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
            if (len(comps) == 4 and
                not self.__current_line().startswith('****')):
                result.append(comps[3][1:])
            self.__current_index += 1
        return result

    def __read_pssms(self):
        """reads and returns the PSSM's"""
        result = []
        self.__find_line_starting_with(' *** Interesting motifs (highest-ranking)')
        self.__next_nonempty_line()
        while (not self.__input_end_reached() and
               not re.match('.*not highest-ranking.*$', self.__current_line())):
            pssm_name = self.__base_name + '_' + self.__current_line().strip()
            result.append(pssm.Pssm(pssm_name, self.__read_frequency_matrix()))
            self.__current_index += 2  # skip over the empty line and the ====
            self.__next_nonempty_line()  # and position to the next PSSM
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
