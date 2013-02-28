# Script to validate a generated output with RegulonDB
# This measures the quality of a cMonkey/Python run by
# comparing the motifs discovered in the run with the
# transcription factors in RegulonDB
import util
import collections
import numpy
import sqlite3
import sys
import argparse

REGULONDB_PSSMS = 'regulondb_pssms.txt'
MEME_PSSMS = 'meme_pssms.txt'

MEME_HEADER = """MEME version 3.0

ALPHABET= ACGT

strands: + -

Background letter frequencies (from dataset with add-one prior applied):
"""

PSSM = collections.namedtuple('PSSM', ['name', 'binding_sites', 'scores'])

def read_pssm(name, text, line_number):
    # 'Total of binding sites: ...'
    line = text[line_number]
    binding_sites = int(line[24:])
    # skip line 'PSSM size: ...' and other sequence stuff
    scores_found = False
    while not scores_found:
        if len(text[line_number]) > 0 and text[line_number][0] == 'a':
            scores_found = True
        else:
            line_number += 1
    a = [float(num) for num in text[line_number][4:].split()]
    c = [float(num) for num in text[line_number + 1][4:].split()]
    g = [float(num) for num in text[line_number + 2][4:].split()]
    t = [float(num) for num in text[line_number + 3][4:].split()]
    scores = numpy.matrix([a, c, g, t])
    scores = numpy.apply_along_axis(lambda a: (a + 0.01) / sum(a + 0.01), 0, scores).T
    return (line_number + 1, PSSM(name, binding_sites, scores))

def read_pssms():
    print "reading PSSMs..."
    pssm_text = util.read_url_cached('http://regulondb.ccg.unam.mx/data/PSSMSet.txt',
                                     'cache/regulondb_pssms.txt').split('\n')
    num_lines = len(pssm_text)
    line_number = 0
    pssms = []
    while line_number < num_lines:
        line = pssm_text[line_number]
        if line.startswith('Transcription Factor Name: '):
            line_number, pssm = read_pssm(line[27:], pssm_text, line_number + 1)
            pssms.append(pssm)
        line_number += 1
    return pssms

def write_pssm_file(pssms):
    with open(REGULONDB_PSSMS, 'w') as outfile:
        outfile.write(MEME_HEADER)
        outfile.write('A 0.283 C 0.217 G 0.217 T 0.283\n')
        for pssm in pssms:
            outfile.write('\n')
            outfile.write('MOTIF %s\n' % pssm.name)
            outfile.write('BL   MOTIF %s width=0 seqs=0\n' % pssm.name)
            outfile.write('letter-probability matrix: alength= 4 w= %s nsites= 20 E= 1.000e+00\n'
                          % pssm.scores.shape[0])
            for row in pssm.scores:
                outfile.write('%5.3f %5.3f %5.3f %5.3f\n' %
                              (row[0], row[1], row[2], row[3]))


def write_motifs(dbname, iteration, seqtype):
    conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    cursor.execute('select rowid, cluster, motif_num from motif_infos ' +
                   'where iteration = ? and seqtype = ?', [iteration, seqtype])
    motifs = cursor.fetchall()

    with open(MEME_PSSMS, 'w') as outfile:
        outfile.write(MEME_HEADER)
        outfile.write('A 0.283 C 0.217 G 0.217 T 0.283\n')

        for motif_id, cluster, motif_num in motifs:
            name = 'MOT_%d_%d' % (cluster, motif_num)
            outfile.write('\n')
            outfile.write('MOTIF %s\n' % name)
            outfile.write('BL   MOTIF %s width=0 seqs=0\n' % name)
            cursor.execute('select a, c, g, t from motif_pssm_rows where motif_info_id = ?' +
                           ' order by row', [motif_id])
            scores = numpy.array(cursor.fetchall())
            outfile.write(('letter-probability matrix: alength= 4 w= %d nsites= 20 ' +
                          'E= 1.000e+00\n') % scores.shape[0])
            for row in scores:
                outfile.write('%5.3f %5.3f %5.3f %5.3f\n' %
                              (row[0], row[1], row[2], row[3]))
        
    cursor.close()
    conn.close()

if __name__ == '__main__':
    description = "tomtom_verify - TomTom verification against RegulonDB"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--dbfile', required=True, help='E.coli cMonkey database')
    args = parser.parse_args()

    pssms = read_pssms()
    print '# RegulonDB PSSMs: ', len(pssms)
    write_pssm_file(pssms)
    print 'writing cMonkey PSSMs'
    write_motifs(args.dbfile, 2001, 'upstream')
