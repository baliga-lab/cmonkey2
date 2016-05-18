# vi: sw=4 ts=4 et:
"""seqtools.py - utilities to operate on genomic sequences

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import itertools as it
import re
import logging
import random
import string
import collections
from cmonkey.util import DelimitedFile

try:
    xrange
except NameError:
    xrange = range


logger = logging.getLogger('seqtools')
logger.setLevel(logging.DEBUG)

# Location and Feature are implemented as named tuples to optimize
# memory and performance
# representation of a genomic position, a simple value object
Location = collections.namedtuple('Location',
                                  ['contig', 'start', 'end', 'reverse'])

# representation of a feature. Just a value object
Feature = collections.namedtuple('Feature',
                                 ['id', 'ftype', 'name', 'location'])


def read_feature(line):
    """Creates and adds a feature and associated contig from current
    DelimitedFile line"""
    contig = line[3]
    # contig is sometimes badly formed with ':'s instead of '_'s
    contig = contig.replace(':', '_')

    is_reverse = False
    if line[6] == 'R':
        is_reverse = True

    # note that feature positions can sometimes start with a '>'
    # or '<', so make sure it is stripped away
    return Feature(line[0], line[1], line[2],
                   Location(contig,
                            int(string.lstrip(line[4], '<>')),
                            int(string.lstrip(line[5], '<>')),
                            is_reverse))


def read_features_from_file(filename):
    """Returns a list containing the features"""
    features = {}
    dfile = DelimitedFile.read(filename, comment='--')
    for line in dfile.lines():
        features[line[0]] = read_feature(line)
    return features


def extract_upstream(source, location, distance):
    """Extract a subsequence of the specified  size from the source sequence
    Depending on the strand orientation, the sequence is cut around either
    the start or the end position"""
    if location.reverse:
        winstart = location.end + 1 + distance[0]
        winend = location.end + 1 + distance[1]
    else:
        winstart = location.start - 1 - distance[1]
        winend = location.start - 1 - distance[0]

    final_location = Location(location.contig, winstart, winend,
                              location.reverse)
    return (final_location,
            subsequence(source, winstart, winend, location.reverse))


def extract_downstream(source, location, distance):
    """Extract a subsequence of the specified  size from the source sequence
    Depending on the strand orientation, the sequence is cut around either
    the start or the end position. NOTE that HERE, distance =
      (number of bases upstream of end, number of bases downstream of end).
    This is the most sensible, but it is NOT the same format as extract_upstream."""
    if location.reverse:
        winstart = location.start + 1 - distance[1]
        winend = location.start + 1 + distance[0]
    else:
        winstart = location.end - 1 - distance[0]
        winend = location.end - 1 + distance[1]

    final_location = Location(location.contig, winstart, winend,
                              location.reverse)
    return (final_location,
            subsequence(source, winstart, winend, location.reverse))


def subsequence(sequence, start, stop, reverse=False):
    """extracts a subsequence from a longer genomic sequence by coordinates.
    If reverse is True, the result string's reverse complement is
    calculated. Not that the start/stop positions are shifted to comply with
    the original cMonkey's behavior
    """
    if start < 1:
        start = 1
    lseq = len(sequence)
    if stop > lseq:
        stop = lseq + 1
    result = sequence[start - 1:stop - 1]
    if reverse:
        result = revcomp(result)
    return result


REV_DICT = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}


def revcomp(sequence):
    """compute the reverse complement of the input string"""
    return "".join([__revchar(c) for c in sequence[::-1]])


def __revchar(nucleotide):
    """for a nucleotide character, return its complement"""
    nucleotide = nucleotide.upper()
    if nucleotide in REV_DICT:
        return REV_DICT[nucleotide]
    else:
        return nucleotide


def subseq_counts(seqs, subseq_len):
    """return a dictionary containing for each subsequence of length
    subseq_len their respective count in the input sequences"""
    counts = {}
    for seq in seqs:
        for index in xrange(0, len(seq) - subseq_len + 1):
            subseq = seq[index:index + subseq_len]
            if not subseq in counts:
                counts[subseq] = 0
            counts[subseq] += 1
    return counts


def subseq_frequencies(seqs, subseq_len):
    """return a dictionary containing for each subsequence of
    length subseq_len their respective frequency within the
    input sequences"""
    result = {}
    counts = subseq_counts(seqs, subseq_len)
    total = sum([count for count in counts.values()])
    for subseq, count in counts.items():
        result[subseq] = float(count) / float(total)
    return result


def markov_background(seqs, order):
    """computes the markov background model of the specified
    order for the given input sequences. This is implemented
    by gathering the frequencies of subsequences of length
    1,..,(order + 1)"""
    result = []
    seqs = replace_degenerate_residues(seqs)
    for subseq_len in xrange(1, (order + 2)):
        result.append(subseq_frequencies(seqs, subseq_len))
    return result


def all_kmers(length, seqs, seq=[], pos=0, choices=['A', 'C', 'G', 'T']):
    """works with lists intead of strings as this is believed to be faster"""
    #logger.debug(pos)
    if seq == []:
        for i in range(length):
            seq.append('')
    for choice in choices:
        seq[pos] = choice
        if pos == length - 1:
            kmer = string.join(seq, '')
            #logger.debug('added kmer %s' %kmer)
            seqs.append(kmer)
        else:
            all_kmers(length, seqs, seq, pos + 1, choices)


def replace_degenerate_residues(seqs):
    """gets rid of funny characters in gene sequences by employing a
    replacement strategy"""
    replacements = {'R': ['G', 'A'], 'Y': ['T', 'C'], 'K': ['G', 'T'],
                    'M': ['A', 'C'], 'S': ['G', 'C'], 'W': ['A', 'T'],
                    'N': ['G', 'A', 'T', 'C'],
                    ' ': [' ']}

    pat = re.compile('[ACGTX]*([^ACGTX])[ACGTX]*')
    result = []
    for seq in seqs:
        seq = seq.strip()  # For some reasons, there were cases with newlines in the beginning
        for match in pat.finditer(seq):
            replace_chars = replacements[seq[match.start(1)]]
            replace_char = replace_chars[random.randint(
                0, len(replace_chars) - 1)]
            seq = seq[:match.start(1)] + replace_char + seq[match.end(1):]
        result.append(seq)
    return result


def read_sequences_from_fasta_string(fasta_string):
    """reads the sequences contained in a FASTA string"""
    lines = fasta_string.split('\n')
    sequences = []
    seqbuffer = ""
    seqname = None
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if len(seqbuffer) > 0:
                sequences.append((seqname, seqbuffer))
                seqbuffer = ""
            seqname = line[1:]
        elif line and len(line) > 0:
            seqbuffer += line
    # add the last line
    if len(seqbuffer) > 0:
        sequences.append((seqname, seqbuffer))
    return sequences


def read_sequences_from_fasta_file(filepath):
    """Read the sequences from the specified FASTA file"""
    with open(filepath) as inputfile:
        fasta_string = inputfile.read()
    return read_sequences_from_fasta_string(fasta_string)


def write_sequences_to_fasta_file(outputfile, seqs):
    """Write a list of sequence tuples to the specified outputfile"""
    for seq in seqs:
        if len(seq[1]) == 0:
            continue  # do not allow 0 length sequences
        outputfile.write('>%s\n' % seq[0])
        outputfile.write('%s\n' % seq[1])


__all__ = ['subsequence', 'extract_upstream', 'markov_background',
           'read_sequences_from_fasta_string',
           'read_sequences_from_fasta_file',
           'write_sequences_to_fasta_file', 'Feature', 'read_features_from_file']
