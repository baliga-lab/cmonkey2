# vi: sw=4 ts=4 et:
"""meme.py - MEME format readers
Handles the file format for MEME for versions >= 4.3.0

This file is part of cm2meme. Please see README and LICENSE for
more information and licensing details.
"""
import subprocess
import logging
import os
import shutil
import re
import collections
from pkg_resources import Requirement, resource_filename, DistributionNotFound

import cmonkey.meme.util as util

try:
    xrange
except NameError:
    xrange = range


class MemeMotifInfo:
    """Only a motif's info line, the
    probability matrix and the site information is relevant"""
    # pylint: disable-msg=R0913
    def __init__(self, pssm, motif_num, width, num_sites, llr, evalue, sites):
        """Creates a MemeMotifInfo instance"""
        self.pssm = pssm
        self.motif_num = motif_num
        self.width = width
        self.num_sites = num_sites
        self.llr = llr
        self.evalue = evalue
        self.sites = sites

    def consensus_string(self, cutoff1=0.7, cutoff2=0.4):
        """returns the consensus string from the pssm table
        remember: letter order is ACGT"""
        alphabet = 'ACGT'
        result = ""
        for row in xrange(len(self.pssm)):
            rowvals = self.pssm[row]
            max_index = rowvals.index(max(rowvals))
            score = rowvals[max_index]
            if score < cutoff2:
                result += 'n'
            elif score < cutoff1:
                result += alphabet[max_index].lower()
            else:
                result += alphabet[max_index]
        return result

    def __repr__(self):
        """returns the string representation"""
        return ("Motif width: %s sites: %s llr: %s e-value: %s" %
                (str(self.width), str(self.num_sites), str(self.llr),
                 str(self.evalue)))


def from_text(output_text, num_motifs):
    """Reads a string in meme output format into a list of
    MemeMotifInfo objects"""

    def extract_width(infoline):
        """extract the width value from the info line"""
        return int(util.extract_regex('width =\s+\d+', infoline))

    def extract_num_sites(infoline):
        """extract the sites value from the info line"""
        return int(util.extract_regex('sites =\s+\d+', infoline))

    def extract_llr(infoline):
        """extract the llr value from the info line"""
        return int(util.extract_regex('llr =\s+\d+', infoline))

    def extract_evalue(infoline):
        """extract the e-value from the info line"""
        return float(util.extract_regex('E-value =\s+\S+', infoline))

    def next_info_line(motif_number, lines, new_motif_names):
        """finds the index of the next info line for the specified motif number
        1-based """
        if new_motif_names:
            return util.next_regex_index('MOTIF [^ ]+ MEME-' + str(motif_number) + '.*',
                                         0, lines)
        else:
            return util.next_regex_index('MOTIF\s+' + str(motif_number) + '.*',
                                         0, lines)

    def next_sites_index(start_index, lines, new_motif_names):
        """returns the next sites index"""
        if new_motif_names:
            return util.next_regex_index('[\t]Motif [^ ]+ MEME-\d+ sites sorted by position ' +
                                         'p-value', start_index, lines)
        else:
            return util.next_regex_index('[\t]Motif \d+ sites sorted by position ' +
                                         'p-value', start_index, lines)

    def read_sites(start_index, lines, new_motif_names):
        """reads the sites"""
        sites_index = next_sites_index(start_index, lines, new_motif_names)
        pattern = re.compile(
            "(\S+)\s+([+-])\s+(\d+)\s+(\S+)\s+(\S+) (\S+) (\S+)?")
        current_index = sites_index + 4
        line = lines[current_index]
        sites = []
        while not line.startswith('----------------------'):
            match = pattern.match(line)
            if match is None:
                logging.error("ERROR in read_sites(), line(#%d) is: '%s'", current_index, line)
            sites.append((match.group(1), match.group(2), int(match.group(3)),
                          float(match.group(4)),
                          match.group(5), match.group(6), match.group(7)))
            current_index += 1
            line = lines[current_index]
        return sites

    def read_pssm(start_index, lines, new_motif_names):
        """reads the PSSM, in this case it's what is called the probability
        matrix in the meme output"""
        pssm_index = next_pssm_index(start_index, lines, new_motif_names)
        current_index = pssm_index + 3
        line = lines[current_index]
        pattern = re.compile("\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")
        rows = []
        while not line.startswith('----------------------'):
            match = pattern.match(line)
            if match is None:
                logging.error("ERROR in read_pssm(), line(#%d) is: '%s'", current_index, line)
            rows.append([float(match.group(1)), float(match.group(2)),
                         float(match.group(3)), float(match.group(4))])
            current_index += 1
            line = lines[current_index]
        return rows

    def next_pssm_index(start_index, lines, new_motif_names):
        """determines the next PSSM start index"""
        if new_motif_names:
            return util.next_regex_index('[\t]Motif [^ ]+ MEME-\d+ position-specific ' +
                                         'probability matrix', start_index, lines)
        else:
            return util.next_regex_index('[\t]Motif \d+ position-specific ' +
                                         'probability matrix', start_index, lines)

    def read_motif_info(motif_number, lines, new_motif_names):
        """Reads the MemeMotifInfo with the specified number from the input"""
        info_line_index = next_info_line(motif_number, lines, new_motif_names)
        info_line = lines[info_line_index]
        return MemeMotifInfo(read_pssm(info_line_index + 1, lines, new_motif_names),
                             motif_number,
                             extract_width(info_line),
                             extract_num_sites(info_line),
                             extract_llr(info_line),
                             extract_evalue(info_line),
                             read_sites(info_line_index + 1, lines, new_motif_names))

    def read_version(lines):
        """Retrieve the file format version.
        Returns: version as a tuple (major, minor, patch)"""
        for line in lines:
            if line.startswith('MEME version'):
                chopped = line.replace('MEME version ', '')
                try:
                    space_pos = chopped.index(' ')
                    return list(map(int, chopped[:space_pos].split('.')))
                except:
                    return list(map(int, chopped.split('.')))
        return None

    def meme_version_conventions(version):
        """Determine conventions used in a file format version
        For now, this only applies to motif naming
        """
        major, minor, patch = version
        new_motif_names = False
        if major == 4:
            if minor == 11:
                new_motif_names = patch >= 4
            elif minor > 11:
                new_motif_names = True
        elif major > 4:
            new_motif_names = True
        return new_motif_names

    lines = output_text.split('\n')
    new_motif_names = meme_version_conventions(read_version(lines))

    result = []
    for motif_number in xrange(1, num_motifs + 1):
        result.append(read_motif_info(motif_number, lines, new_motif_names))
    return result



######################################################################
### Export
###############################3

MEME_FILE_HEADER = """MEME version 3.0

ALPHABET= ACGT

strands: + -

Background letter frequencies (from dataset with add-one prior applied):
A %.3f C %.3f G %.3f T %.3f
"""

"""
def write_pssm(session, outfile, motif_info_id, evalue, num_sites):
    # writes a single PSSM to the given file
    outfile.write('\nMOTIF %d\n' % motif_info_id)
    outfile.write('BL   MOTIF %s width=0 seqs=0\n' % motif_info_id)

    pssm_rows = [(row.a, row.c, row.g, row.t)
        for row in session.query(cm2db.MotifPSSMRow).filter(cm2db.MotifPSSMRow.motif_info_id == motif_info_id)]

    outfile.write('letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e\n' % (len(pssm_rows), num_sites, evalue))
    for a, c, g, t in pssm_rows:
        outfile.write('%5.3f %5.3f %5.3f %5.3f\n' % (a, c, g, t))


def write_motifs2meme(session, filepath):
    # Write the motifs to a MEME file and returns True if successful
    #Currently, this only works if there is global background data in the database
    freqs = {gb_row.subsequence: gb_row.pvalue
        for gb_row in session.query(cm2db.GlobalBackground).filter(cm2db.GlobalBackground.subsequence.in_(['A', 'C', 'G', 'T']))}

    if len(freqs) >= 4 and 'A' in freqs and 'C' in freqs and 'G' in freqs and 'T' in freqs:
        logging.debug('retrieving letter frequency from global background distribution')
        with open(filepath, 'w') as outfile:
            outfile.write(MEME_FILE_HEADER % (freqs['A'], freqs['C'], freqs['G'], freqs['T']))

            iteration = session.query(func.max(cm2db.MotifInfo.iteration))
            for mi in session.query(cm2db.MotifInfo).filter(cm2db.MotifInfo.iteration == iteration):
                if mi.num_sites > 0:
                    write_pssm(session, outfile, mi.rowid, mi.evalue, mi.num_sites)
                else:
                    # if we don't run MEME, but Weeder, num_sites is 0
                    write_pssm(session, outfile, mi.rowid, mi.evalue, mi.num_annotations)

        return True
    else:
        logging.warn('no global background distribution found')
        return False


EVALUE_CUTOFF = 100
RESID_CUTOFF  = 0.8
DIST_METHOD   = "ed"
Q_THRESHOLD   = 0.5
MIN_OVERLAP   = 4
MOTIF_PSEUDO  = 0.0


def run_tomtom(session, targetdir, version, q_thresh=Q_THRESHOLD, dist_method=DIST_METHOD,
               min_overlap=MIN_OVERLAP, motif_pseudo=MOTIF_PSEUDO):
    # a wrapper around the tomtom script
    targetfile = os.path.join(targetdir, 'post.tomtom.meme')
    queryfile = targetfile
    if write_motifs2meme(session, targetfile):
        if version.startswith('4.11'):
            # tomtom command parameter names changed in version 4.11.x:
            #   - target-pseudo and query-pseudo merged into motif-pseudo
            #   - q-thresh changed to thresh
            command = ['tomtom',
                       '-verbosity', '1',
                       '-thresh', '%.3f' % q_thresh,
                       '-dist', dist_method,
                       '-min-overlap', '%d' % min_overlap,
                       '-text',
                       '-motif-pseudo', '%.3f' % motif_pseudo]
        else:
            command = ['tomtom',
                       '-verbosity', '1',
                       '-q-thresh', '%.3f' % q_thresh,
                       '-dist', dist_method,
                       '-min-overlap', '%d' % min_overlap,
                       '-text',
                       '-query-pseudo', '%.3f' % motif_pseudo,
                       '-target-pseudo', '%.3f' % motif_pseudo]

        logging.debug(" ".join(command))

        # Tomtom versions > 4.8.x drop the target and query switches
        if version == '4.3.0':
            command.extend(['-target', targetfile, '-query', queryfile])
        else:
            command.extend([targetfile, queryfile])

        try:
            output = subprocess.check_output(command).decode('utf-8')
            lines = output.split('\n')[1:]
            for line in lines:
                if len(line.strip()) > 0:
                    row = line.strip().split('\t')
                    motif1 = int(row[0])
                    motif2 = int(row[1])
                    if motif1 != motif2:
                        pvalue = float(row[3])
                        session.add(cm2db.TomtomResult(motif_info_id1=motif1, motif_info_id2=motif2, pvalue=pvalue))
            session.commit()
        except:
            raise
"""

__all__ = ['read_meme_output']
