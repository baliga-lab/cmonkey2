# vi: sw=4 ts=4 et:
"""pssm.py - cMonkey pssm module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import math


TWO_BASE_LETTERS = ['Y', 'R', 'W', 'S', 'K', 'M']
THREE_BASE_LETTERS = ['V', 'H', 'D', 'B']


class Pssm(object):
    """A PSSM class to interface with Weeder"""
    __slots__ = ('name', 'values', 'sites', 'e_value')

    def __init__(self, name, values=None, e_value=None, sites=None):
        """Create a PSSM object"""
        self.name = name
        if values:
            self.values = values
        else:
            self.values = []
        self.sites = sites
        self.e_value = e_value

    def __getitem__(self, row):
        """accesses the row at the specified index"""
        return self.values[row]

    def sequence_length(self):
        """returns the sequence length this PSSM is based on"""
        return len(self.values)

    def to_logodds_string(self, at_freq=0.25, cg_freq=0.25):
        """returns a motif representation in log-odds format"""
        def log_odds(pvalue, freq):
            """returns the log-odds value"""
            if pvalue == 0.0:
                return int(round(math.log(1e-2 / freq, 2)))  # was 1e-300
            else:
                return int(round(math.log(pvalue / freq, 2)))

        result = ('log-odds matrix: alength= 4 w= %d\n' %
                  self.sequence_length())
        for row in self.values:
            result += ('%7d %7d %7d %7d\n' %
                       (log_odds(row[0], at_freq),
                        log_odds(row[1], cg_freq),
                        log_odds(row[2], cg_freq),
                        log_odds(row[3], at_freq)))
        return result

    def consensus_motif(self, limit1=0.6, limit2=0.8, three=True):
        """returns the consensus motif"""
        def column_consensus(row):
            """returns the column consensus"""
            if self.values[row][0] >= limit1:
                return 'A'
            elif self.values[row][1] >= limit1:
                return 'C'
            elif self.values[row][2] >= limit1:
                return 'G'
            elif self.values[row][3] >= limit1:
                return 'T'
            else:
                return compute_column_consensus(row)

        def compute_column_consensus(row):
            """computes column consensus"""
            two_base = [
                self.values[row][1] + self.values[row][3],
                self.values[row][0] + self.values[row][2],
                self.values[row][0] + self.values[row][3],
                self.values[row][1] + self.values[row][2],
                self.values[row][2] + self.values[row][3],
                self.values[row][0] + self.values[row][1]]

            three_base = [
                self.values[row][0] + self.values[row][1] +
                self.values[row][2],
                self.values[row][0] + self.values[row][1] +
                self.values[row][3],
                self.values[row][0] + self.values[row][2] +
                self.values[row][3],
                self.values[row][1] + self.values[row][2] +
                self.values[row][3]]

            pmax = 0.0
            max_letter = 'N'
            for row in range(len(two_base)):
                if two_base[row] > pmax:
                    pmax = two_base[row]
                    max_letter = TWO_BASE_LETTERS[row]
            if pmax <= limit2 and three:
                for row in range(len(three_base)):
                    if three_base[row] > pmax:
                        pmax = three_base[row]
                        max_letter = THREE_BASE_LETTERS[row]
            if pmax <= limit2:
                max_letter = 'N'
            return max_letter

        result = ""
        for col in range(len(self.values)):
            result += column_consensus(col)
        return result


def read_fasta(infile):
    """creates a list of Pssm objects from a FASTA file"""
    def next_start_index(lines, starting_at):
        """returns the next start index of an entry starting at
        a specific index"""
        current_index = starting_at
        while current_index < len(lines):
            line = lines[current_index]
            if len(line) > 0 and line[0] == '>':
                return current_index
            current_index += 1
        return -1

    def read_values(line):
        """returns the float values in the specified line"""
        return [float(strvalue) for strvalue in line.strip().split(' ')]

    def read_pssm(lines, starting_at):
        """reads a PSSM from the current line index"""
        elem_index = next_start_index(lines, starting_at)
        if elem_index >= 0:
            line = lines[elem_index]
            name = line[1:].strip()
            matrix = []
            for offset in range(1, 5):
                matrix.append(read_values(lines[elem_index + offset]))

            values = [[matrix[col][row] for col in range(4)]
                      for row in range(len(matrix[0]))]
            return Pssm(name, values), elem_index + 5
        else:
            return None, -1

    lines = infile.readlines()
    next_index = 0
    result = []
    while next_index >= 0 and next_index < len(lines):
        pssm, next_index = read_pssm(lines, next_index)
        if pssm is not None:
            result.append(pssm)
    return result
