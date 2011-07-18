"""thesaurus.py - cMonkey thesaurus module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


class DelimitedFileThesaurusFactory:
    """creates a thesaurus from a delimited file"""

    def __init__(self, dfile):
        """creates an instance of this factory"""
        self.dfile = dfile

    def create(self):
        """returns a thesaurus dictionary"""
        result = {}
        for line in self.dfile.lines():
            result[line[0]] = line[1]
        return result
