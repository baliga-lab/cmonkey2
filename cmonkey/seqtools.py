"""seqtools.py - utilities to operate on genomic sequences

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


class Location:  # pylint: disable-msg=R0903
    """representation of a genomic position, a simple value object"""
    def __init__(self, contig, start, end, reverse):
        """Creates a location instance"""
        self.contig = contig
        self.start = start
        self.end = end
        self.reverse = reverse

    def __repr__(self):
        """String representation"""
        return ("contig: %s s: %d e: %d rev: %s" %
                (self.contig, self.start,
                 self.end, str(self.reverse)))

    def __eq__(self, other):
        """implements the == operation"""
        return (self.contig == other.contig and self.start == other.start and
                self.end == other.end and self.reverse == other.reverse)

    def __ne__(self, other):
        """implements the != operation"""
        return not self.__eq__(other)


class Feature:  # pylint: disable-msg=R0902
    """representation of a feature. Just a value object"""

    def __init__(self, feature_id, feature_type, name, location):
        """Create a Feature instance"""
        # pylint: disable-msg=R0913
        self.__feature_id = feature_id
        self.__feature_type = feature_type
        self.__name = name
        self.__location = location

    def id(self):  # pylint: disable-msg=C0103
        """returns the feature id"""
        return self.__feature_id

    def location(self):
        """returns this feature's location"""
        return self.__location

    def __repr__(self):
        """returns the string representation"""
        return ("%s[%s] - %s, %s" %
                (self.__feature_id, self.__feature_type,
                 self.__name, repr(self.__location)))


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


def subsequence(sequence, start, stop, reverse=False):
    """extracts a subsequence from a longer genomic sequence by coordinates.
    If reverse is True, the result string's reverse complement is
    calculated. Not that the start/stop positions are shifted to comply with
    the original cMonkey's behavior
    """
    result = sequence[start - 1:stop - 1]
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


def subseq_counts(seqs, subseq_len):
    """return a dictionary containing for each subsequence of length
    subseq_len their respective count in the input sequences"""
    counts = {}
    for seq in seqs:
        for index in range(0, len(seq) - subseq_len + 1):
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
    for subseq_len in range(1, (order + 2)):
        result.append(subseq_frequencies(seqs, subseq_len))
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
        outputfile.write('>%s\n' % seq[0])
        outputfile.write('%s\n' % seq[1])

__all__ = ['subsequence', 'extract_upstream', 'markov_background',
           'read_sequences_from_fasta_string',
           'read_sequences_from_fasta_file',
           'write_sequences_to_fasta_file', 'Feature']
