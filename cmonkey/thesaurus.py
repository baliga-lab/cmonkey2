# vi: sw=4 ts=4 et:
"""thesaurus.py - cMonkey thesaurus module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import re
import cmonkey.util as util


def create_from_delimited_file1(dfile):
    """creates a thesaurus from a delimited file where the format is
    <alternative>SEPARATOR<original>
    ..."""
    return {line[0]: line[1] for line in dfile.lines}


def create_from_delimited_file2(dfile, case_sensitive):
    """creates a thesaurus from a delimited file where the format is
    <original>SEPARATOR<alt1>;<alt2>;...
    ..."""
    def fix_case(s):
        return s if case_sensitive else s.upper()

    if isinstance(dfile, str):
        dfile = util.read_dfile(dfile, sep=',', has_header=False)
    result = {}
    for line in dfile.lines:
        original = fix_case(line[0])  # original should map to itself
        result[original] = original
        for alternative in line[1].split(';'):
            result[fix_case(alternative)] = original
    return result


def create_from_rsat_feature_names(dfile, key_transforms=None):
    """Uses a feature_names.tab file from RSAT to create a Thesaurus,
    represented as a dictionary.
    key_transforms is an optional list of key modifiers that are applied in
    sequence to the key before it is added to the dictionary
    a key_transform is a function from string -> [string], so multiple
    results can be returned from the transformer
    e.g. the feature_names.tab file stores VNG names with a modification
    suffix that can be removed
    """
    result = {}
    for line in dfile.lines:
        key = line[1]
        alternative = line[0]
        if key_transforms:
            for transform in key_transforms:
                for transform_key in transform(key):
                    result[transform_key] = alternative
        else:
            result[key] = alternative
    return result


def strip_vng_modification(gene):
    """strips 'm' modifier off a VNG name"""
    if re.match('VNG\d{4}.m$', gene):
        return [gene, gene.rstrip('m')]
    else:
        return [gene]


__all__ = ['create_from_delimited_file1', 'create_from_delimited_file2',
           'create_from_rsat_feature_names']
