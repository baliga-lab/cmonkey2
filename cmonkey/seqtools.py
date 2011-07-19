"""seqtools.py - utilities to operate on genomic sequences

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


def extract_upstream(source, start, end, reverse, distance):
    """rule for upstream extraction"""
    return subsequence(source, start + (distance[0] - 1),
                       end + (distance[1] - 1), reverse)


def subsequence(sequence, start, stop, reverse=False):
    """extracts a subsequence from a longer genomic sequence by coordinates.
    If reverse is True, the result string's reverse complement is
    calculated
    """
    result = sequence[start:(stop + 1)]
    if reverse:
        result = revcomp(result)
    return result


def revcomp(sequence):
    """compute the reverse complement of the input string"""
    return "".join([revchar(c) for c in sequence[::-1]])


def revchar(nucleotide):
    """for a nucleotide character, return its complement"""
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'T':
        return 'A'
    else:
        raise ValueError('unknown param: %s' % str(nucleotide))


__all__ = ['subsequence', 'extract_upstream']
