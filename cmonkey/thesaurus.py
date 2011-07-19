"""thesaurus.py - cMonkey thesaurus module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


def create_from_delimited_file1(dfile):
    """creates a thesaurus from a delimited file where the format is
    <alternative>SEPARATOR<original>
    ..."""
    result = {}
    for line in dfile.lines():
        result[line[0]] = line[1]
    return result


def create_from_delimited_file2(dfile):
    """creates a thesaurus from a delimited file where the format is
    <original>SEPARATOR<alt1>;<alt2>;...
    ..."""
    result = {}
    for line in dfile.lines():
        for alternative in line[1].split(';'):
            result[alternative] = line[0]
    return result


def create_from_rsat_feature_names(dfile):
    """Uses a feature_names.tab file from RSAT to create a Thesaurus"""
    result = {}
    for line in dfile.lines():
        result[line[1]] = line[0]
    return result


__all__ = ['create_from_delimited_file1', 'create_from_delimited_file2',
           'create_from_rsat_feature_names']
