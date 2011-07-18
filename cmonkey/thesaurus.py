"""thesaurus.py - cMonkey thesaurus module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


class DelimitedFileThesaurusFactory1:
    """creates a thesaurus from a delimited file where the format is
    <alternative>SEPARATOR<original>
    ...
    """

    def __init__(self, dfile):
        """creates an instance of this factory"""
        self.dfile = dfile

    def create(self):
        """returns a thesaurus dictionary"""
        result = {}
        for line in self.dfile.lines():
            result[line[0]] = line[1]
        return result


class DelimitedFileThesaurusFactory2:
    """creates a thesaurus from a delimited file where the format is
    <original>SEPARATOR<alt1>;<alt2>;...
    ...
    """

    def __init__(self, dfile):
        """creates an instance of this factory"""
        self.dfile = dfile

    def create(self):
        """returns a thesaurus dictionary"""
        result = {}
        for line in self.dfile.lines():
            for alternative in line[1].split(';'):
                result[alternative] = line[0]
        return result
