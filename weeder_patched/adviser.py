#!/opt/local/bin/python
import sys
import seqtools as st
import re
import math

ALPHABET           = 'ACGT'
MAX_PATTERN_LENGTH = 15
NUM_BUCKETS        = 4
ALPHABET_SIZE      = 4

def log(value):
    """This is a wrapper around math.log() to return negative infinity if
    value is 0
    """
    if value == 0:
        return float('-inf')
    else:
        return math.log(value)

class MixPattern:
    def __init__(self, rank, pattern, score, error):
        # base information from mix file
        self.rank       = int(rank)
        self.pattern    = pattern
        self.score      = float(score)
        self.error      = int(error)

        # additional metrics computed
        self.hor        = 0
        self.vert       = 0
        self.hamming    = 0
        self.min        = 0.0
        self.max        = 0.0
        self.counter    = 0
        self.brothers   = []

        self.overlap    = False
        self.inside     = False
        self.outside    = False
        self.matrix = [[0, 0, 0, 0] for i in range(MAX_PATTERN_LENGTH)]

    def usable(self):
        return ((self.vert > 0 and (self.hor > 0 or self.vert > 2)) or
                (self.rank == 1 and len(self.pattern) > 6))

    def __str__(self):
        return "%d. '%s' score = %f, error = %d" % (self.rank,
                                                    self.pattern,
                                                    self.score,
                                                    self.error)
    def __repr__(self):
        return self.__str__()


def read_mixfile(filename):
    pattern = re.compile('^(\d+)\) ([AGTC]+) (\d+\.\d+) (\d+)$')
    with open(filename) as infile:
        lines = infile.readlines()
    entries = []
    for line in lines:
        matcher = pattern.match(line.strip())
        if matcher != None:
            entries.append(MixPattern(matcher.group(1),
                                      matcher.group(2),
                                      matcher.group(3),
                                      matcher.group(4)))
    return entries


def add_to_bucket(bucket, new_entry):
    found = False
    for entry in bucket:
        if (entry.pattern == new_entry.pattern and
            (entry.rank == new_entry.rank and entry.score < new_entry.score or
             entry.rank < new_entry.rank)):
            entry.score = new_entry.score
            entry.error = new_entry.error
            found = True
            break
    if not found:
        bucket.append(new_entry)


def make_buckets(mixentries):
    buckets = [[] for i in range(NUM_BUCKETS)]
    for entry in mixentries:
        bucket_num = (len(entry.pattern) - 6) / 2
        add_to_bucket(buckets[bucket_num], entry)
    return buckets

def compute_pattern_relationships(buckets, reverse):
    for j in range(len(buckets)):
        for tmppat in buckets[j]:
            for i in range(j, len(buckets)):
                for seekpat in buckets[i]:
                    if tmppat.pattern != seekpat.pattern:
                        overlap = st.overlap(seekpat.pattern, tmppat.pattern, reverse)
                        if overlap and tmppat.rank < seekpat.rank:
                            tmppat.vert         = 1
                            tmppat.overlap      = True
                            tmppat.brothers.append(seekpat.pattern)

                            seekpat.vert        = 1
                            seekpat.overlap     = True
                            seekpat.brothers.append(tmppat.pattern)

                        if (st.hamming_distance(seekpat.pattern, tmppat.pattern, reverse) <=
                            seekpat.error) and tmppat.rank < seekpat.rank:
                            tmppat.hamming     = 1
                            tmppat.vert        = 1
                            tmppat.brothers.append(seekpat.pattern)

                            seekpat.hamming     = 1
                            seekpat.vert        = 1
                            seekpat.brothers.append(tmppat.pattern)

                        if (len(seekpat.pattern) != len(tmppat.pattern) and
                            st.inside(tmppat.pattern, seekpat.pattern, reverse)):
                            tmppat.hor         = 1
                            tmppat.inside      = True
                            tmppat.brothers.append(seekpat.pattern)

                            seekpat.hor         = 1
                            seekpat.outside     = True
                            seekpat.brothers.append(tmppat.pattern)


def search_patterns(buckets, seqs, reverse):
    for i in range(len(buckets)):
        for seekpat in buckets[i]:
            maxerr    = seekpat.error
            howmany   = 0
            lastinseq = -1
            expected  = 0

            if seekpat.usable():
                for uu in range(len(seqs)):
                    bestpos = -1
                    # TODO: reverse ignored for now (maah is different in that case)
                    maah = len(seqs[uu][1]) - len(seekpat.pattern) + 1

                    for j in range(maah):
                        k = 0
                        flagnope = False
                        nuclei    = []

                        for pp in range(len(seekpat.pattern)):
                            if st.char_to_int(seqs[uu][1][j + pp]) < 0:
                                flagnope = True

                        if not flagnope:
                            for pp in range(len(seekpat.pattern)):
                                nuclei.append(st.char_to_int(seqs[uu][1][j + pp]))
                                if seqs[uu][1][j + pp].lower() != seekpat.pattern[pp].lower():
                                    k += 1
                        else:
                            k = maxerr + 1
                        
                        # TODO: ignore reverse for now (condition changes in that case)
                        if k <= maxerr:
                            howmany += 1
                            nn = len(seekpat.pattern) - 1
                            for jj in range(nn + 1):
                                seekpat.matrix[jj][nuclei[jj]] += 1

                            if lastinseq != uu:
                                lastinseq = uu
                                bestpos = j
                            k = 0

                seekpat.min = 0
                seekpat.max = 0
                
                for pp in range(len(seekpat.pattern)):
                    maxtmp = 0.0
                    mintmp = 1.0
                    for qq in range(ALPHABET_SIZE):
                        pval = float(seekpat.matrix[pp][qq] / float(howmany))
                        if pval > maxtmp:
                            maxtmp = pval
                        if pval < mintmp:
                            mintmp = pval
                    seekpat.min += log(mintmp + 0.001)
                    seekpat.max += log(maxtmp)
                seekpat.counter = howmany
            """
            print "SEEKPAT '%s' MIN = %f MAX = %f COUNT = %d" % (seekpat.pattern, seekpat.min,
                                                                 seekpat.max, seekpat.counter)
            for pp in range(len(seekpat.pattern)):
                print "[",
                for qq in range(ALPHABET_SIZE):
                    print seekpat.matrix[pp][qq],
                print "]"
                """


def output_pattern(htmlfile, seekpat, seqs, reverse):
    htmlfile.write('<div>%s</div>' % seekpat.pattern);
    if reverse:
        htmlfile.write('<div>%s</div>' % seekpat.pattern[::-1])
    htmlfile.write('<div>%d redundant motifs found:</div>' % len(seekpat.brothers))
    htmlfile.write('<div>%s</div>' % ' - '.join(seekpat.brothers))

    # Print the occurrences table
    htmlfile.write('<div>Best occurrences (match percentage)</div>')
    htmlfile.write('<table><tr><th>Sequence</th><th>Strand</th><th>Oligo</th>' +
                   '<th>Position</th><th>Match</th><tr>')
    maxerr    = seekpat.error
    newmatrix = [[0, 0, 0, 0] for m in range(MAX_PATTERN_LENGTH)]
    howmany   = 0
    lastinseq = -1
    expected  = 0

    for uu in range(len(seqs)):
        bestpos = -1
        # TODO: reverse is ignored for now (maah is dependent on that)
        maah = len(seqs[uu][1]) - len(seekpat.pattern) + 1
        
        for j in range(maah):
            k = 0
            flagnope = False
            nuclei = []

            for pp in range(len(seekpat.pattern)):
                if st.char_to_int(seqs[uu][1][j + pp]) < 0:
                    flagnope = True
            if not flagnope:
                for pp in range(len(seekpat.pattern)):
                    nuclei.append(st.char_to_int(seqs[uu][1][j + pp]))
            else:
                k = maxerr + 1

            # TODO: reverse is ignored for now
            if k <= maxerr:
                nn = len(seekpat.pattern) - 1
                tmpflo = 0.0
                for jj in range(nn + 1):
                    tmpflo += (log(float(seekpat.matrix[jj][nuclei[jj]])) -
                               log(float(seekpat.counter)))

                percentage = (100.0 * (tmpflo - seekpat.min) / (seekpat.max - seekpat.min))
                if percentage > 85:
                    howmany += 1
                    htmlfile.write('<tr>')
                    # TODO: reverse is ignored for now
                    htmlfile.write('<td>%d</td>' % (uu + 1))
                    htmlfile.write('<td>+</td>')
                    htmlfile.write('<td>')
                    if percentage < 90:
                        htmlfile.write('[')
                    else:
                        htmlfile.write('.')

                    for jj in range(nn + 1):
                        htmlfile.write('%c' % ALPHABET[nuclei[jj]])
                        if percentage >= 90:
                            newmatrix[jj][nuclei[jj]] += 1

                    if percentage < 90:
                        htmlfile.write(']')
                    else:
                        htmlfile.write('.')

                    htmlfile.write('</td>')
                    # TODO: reversed is ignored for now
                    htmlfile.write('<td>%d</td>' % (j + 1))

                    htmlfile.write('<td>(%.2f)</td>' % percentage)

                    htmlfile.write('</tr>\n')

    htmlfile.write('</table>')

    # print frequency matrix
    newmin = 0.0
    newmax = 0.0
    for pp in range(len(seekpat.pattern)):
        maxtmp = 0.0
        mintmp = 1.0
        for qq in range(ALPHABET_SIZE):
            pval = float(newmatrix[pp][qq]) / howmany
            if pval > maxtmp:
                maxtmp = pval
            if pval < mintmp:
                mintmp = pval
        newmin += log(mintmp + 0.001)
        newmax += log(maxtmp)
    htmlfile.write("""<table><caption>Frequency Matrix</caption>
<tr><th>&nbsp;</th><th colspan="5">All Occs</th><th>&nbsp;</th><th>Best Occs</th></tr>
<tr><th>&nbsp;</th><th>A</th><th>C</th><th>G</th><th>T</th>
<th>&nbsp;</th><th>A</th><th>C</th><th>G</th><th>T</th></tr>
""")
    for pp in range(len(seekpat.pattern)):
        htmlfile.write('<tr><td>%d</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td>' %
                       (pp + 1, seekpat.matrix[pp][0],
                        seekpat.matrix[pp][1],
                        seekpat.matrix[pp][2],
                        seekpat.matrix[pp][3]))
        htmlfile.write('<td>|</td><td>%d</td><td>%d</td><td>%d</td><td>%d</td></tr>' %
                       (newmatrix[pp][0],
                        newmatrix[pp][1],
                        newmatrix[pp][2],
                        newmatrix[pp][3]))
    htmlfile.write("</table>")

def output_results(seqs, mix_entries, buckets, basename, reverse):
    patlimit   = 4
    has_advice = False

    htmlfilename = '%s.py.html' % basename
    print "Writing output to '%s'" % htmlfilename

    html_head = """<!doctype html><html><head><title>Adviser output</title>
<link rel="stylesheet" href="adviser-style.css">
</head>\n"""
    with open(htmlfilename, 'w') as htmlfile:
        htmlfile.write(html_head)
        if len(mix_entries) == 0:
            print "No motifs"
        else:
            htmlfile.write("<h3>Your sequences:</h3>")

        for i in range(len(seqs)):
            htmlfile.write("Sequence %d : %s<br>" % (i + 1, seqs[i][0]))

        htmlfile.write('<h2>My Advice</h2>')

        if len(buckets[2]) <= 0:
            patlimit = 2
        for rankcount in range(len(buckets[0]) / 2):
            for i in range(4):
                if rankcount == 0 and i == 0:
                    htmlfile.write('<h3>Interesting motifs (highest-ranking) seem to be:</h3>')
                    has_advice = False
                if rankcount == 1 and i == 0:
                    htmlfile.write('<h3>Interesting motifs (not highest-ranking) can also be:</h3>')
                    has_advice = False

                if len(buckets[i]) > 0:
                    pat_index = 0
                    seekpat = buckets[i][pat_index]
                    if seekpat.score > 0:
                        pat_index += 1
                        if seekpat.usable():
                            output_pattern(htmlfile, seekpat, seqs, reverse)
        htmlfile.write("</body></html>")


if __name__ == '__main__':
    print "Adviser Python"
    if len(sys.argv) <= 1:
        print "usage: python adviser.py <fasta-file> [S]"
    else:
        basename = sys.argv[1]
        reverse = False
        print "Processing '%s'" % basename
        seqs = st.read_sequences_from_fasta_file(basename)
        mixfile  = '%s.mix' % basename
        print "Running adviser on '%s'" % mixfile

        mix_entries = read_mixfile(mixfile)  # seems correct
        buckets = make_buckets(mix_entries)  # seems correct
        compute_pattern_relationships(buckets, reverse)  # seems correct
        search_patterns(buckets, seqs, reverse)  # fixed
        output_results(seqs, mix_entries, buckets, sys.argv[1], reverse)
