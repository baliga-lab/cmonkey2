"""pssm.py - cMonkey pssm module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


TWO_BASE_LETTERS   = ['Y', 'R', 'W', 'S', 'K', 'M']
THREE_BASE_LETTERS = ['V', 'H', 'D', 'B']


class Pssm:
    """A PSSM class to interface with Weeder"""

    def __init__(self, name, values=[], num_sites=0, e_value=None):
        """Create a PSSM object"""
        self.__name = name
        self.__values = values
        self.__num_sites = num_sites
        self.__e_value = e_value

    def __getitem__(self, row):
        """accesses the row at the specified index"""
        return self.__values[row]

    def name(self):
        """returns the name"""
        return self.__name

    def num_sites(self):
        """returns the number of sites"""
        return self.__num_sites

    def e_value(self):
        """returns the e-value"""
        return self.__e_value

    def consensus_motif(self, limit1=0.6, limit2=0.8, three=True):
        """returns the consensus motif"""
        def column_consensus(row):
            """returns the column consensus"""
            if self.__values[row][0] >= limit1:
                return 'A'
            elif self.__values[row][1] >= limit1:
                return 'C'
            elif self.__values[row][2] >= limit1:
                return 'G'
            elif self.__values[row][3] >= limit1:
                return 'T'
            else:
                return compute_column_consensus(col)

        def compute_column_consensus(row):
            """computes column consensus"""
            two_base = [
                self.__values[row][1] + self.__values[row][3],
                self.__values[row][0] + self.__values[row][2],
                self.__values[row][0] + self.__values[row][3],
                self.__values[row][1] + self.__values[row][2],
                self.__values[row][2] + self.__values[row][3],
                self.__values[row][0] + self.__values[row][1]]

            three_base = [
                self.__values[row][0] + self.__values[row][1] +
                self.__values[row][2],
                self.__values[row][0] + self.__values[row][1] +
                self.__values[row][3],
                self.__values[row][0] + self.__values[row][2] +
                self.__values[row][3],
                self.__values[row][1] + self.__values[row][2] +
                self.__values[row][3]]

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
                maxletter = 'N'
            return maxletter

        result = ""
        for col in range(len(self.__values)):
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
        if pssm != None:
            result.append(pssm)
    return result
