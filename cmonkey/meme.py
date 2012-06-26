# vi: sw=4 ts=4 et:
"""meme.py - cMonkey meme suite integration
These are functions and classes to conveniently call the needed
commands of the MEME suite in order to find motifs

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import subprocess
import tempfile
import logging
import seqtools as st
import os
import util
import re


class MemeRunResult:
    """Result data for a single MEME run"""
    def __init__(self, pe_values, annotations, motif_infos):
        """constructor"""
        self.pe_values = pe_values
        self.annotations = annotations
        self.motif_infos = motif_infos


class MemeSuite:
    """Regard the meme suite as a unit of tools. This helps
    us capturing things like versions, global settings and data
    passing
    These represent the command line tools and currently are specific
    for 4.3.0, in 4.6.1 the interface has changed, later the MemeSuite
    will take a MemeSuiteVersion object to make these changes transparent.

    dust - remove low-complexity regions or sequence repeats
    meme - discover motifs in a set of sequences
    mast - search for a group of motifs in a set of sequences
    """
    def __init__(self, max_width=24, use_revcomp=True, background_file=None,
            remove_tempfiles=True):
        """Create MemeSuite instance"""
        self.__max_width = max_width
        self.__use_revcomp = use_revcomp
        self.__background_file = background_file
        self.__remove_tempfiles = remove_tempfiles

    def global_background_file(self):
        """returns the global background file used with this meme suite
        instance"""
        return self.__background_file

    def max_width(self):
        """returns the max_width attribute"""
        return self.__max_width

    def remove_low_complexity(self, seqs):
        """send sequences through dust filter, send only those
        to dust that are larger than max_width"""
        def process_with_dust(seqs):
            """data conversion from and to dust tool"""
            dust_tmp_file = None
            with tempfile.NamedTemporaryFile(prefix='dust',
                                             delete=False) as dust_input:
                for feature_id, seq in seqs.items():
                    dust_input.write(">%s\n" % feature_id)
                    dust_input.write("%s\n" % seq[1])
                dust_tmp_file = dust_input.name
                #logging.info("DUST input written to: %s", dust_input.name)
            seqpairs = st.read_sequences_from_fasta_string(
                self.dust(dust_tmp_file))
            os.remove(dust_tmp_file)
            result = {}
            for feature_id, seq in seqpairs:
                result[feature_id] = seq
            return result

        seqs_for_dust = {}
        for feature_id, seq in seqs.items():
            if len(seq[1]) > self.__max_width:
                seqs_for_dust[feature_id] = seq
        # only non-empty-input gets into dust, dust can not
        # handle empty input
        if len(seqs_for_dust) > 0:
            return process_with_dust(seqs_for_dust)
        else:
            return {}

    def __call__(self, params):
        """Runs the meme tool. input_seqs is a dictionary of
        (feature_id : (location, sequence)) that are to be provided as meme
        input, all_seqs is a dictionary that provides all sequences used
        in the cMonkey run, which will be used to compute background
        distribution.
        Note: To more closely resemble the original R algorithm, we provide
        ----- the sorted feature ids so MEME will return the same output"""
        input_seqs = params.seqs
        all_seqs = params.used_seqs
        
        def background_file():
            """decide whether to use global or specific background file"""
            if self.__background_file != None:
                logging.info("using global background")
                return self.__background_file
            else:
                bgseqs = {feature_id: all_seqs[feature_id]
                          for feature_id in all_seqs
                          if feature_id not in params.feature_ids}
                return make_background_file(bgseqs, self.__use_revcomp)

        #logging.info("run_meme() - # seqs = %d", len(input_seqs))
        bgfile = background_file()
        #logging.info("created background file in %s", bgfile)
        seqfile = self.make_sequence_file(
            [(feature_id, input_seqs[feature_id])
             for feature_id in params.feature_ids if feature_id in input_seqs])
        #logging.info("created sequence file in %s", seqfile)
        motif_infos, output = self.meme(seqfile, bgfile, params.num_motifs,
                                        params.previous_motif_infos)

        # run mast
        meme_outfile = None
        with tempfile.NamedTemporaryFile(prefix='meme.out.',
                                         delete=False) as outfile:
            meme_outfile = outfile.name
            outfile.write(output)
        logging.info('wrote meme output to %s', meme_outfile)
        dbfile = self.make_sequence_file(
            [(feature_id, locseq[1])
             for feature_id, locseq in all_seqs.items()])
        logging.info('created mast database in %s', dbfile)
        try:
            mast_output = self.mast(meme_outfile, dbfile, bgfile)
            pe_values, annotations = read_mast_output(mast_output,
                                                      input_seqs.keys())
            return MemeRunResult(pe_values, annotations, motif_infos)
        except:
            return MemeRunResult([], {}, [])
        finally:
            if self.__remove_tempfiles:
                #logging.info("DELETING ALL TMP FILES...")
                try:
                    os.remove(seqfile)
                except:
                    logging.warn("could not remove tmp file: '%s'", seqfile)
                try:
                    os.remove(meme_outfile)
                except:
                    logging.warn("could not remove tmp file: '%s'", meme_outfile)
                try:
                    os.remove(dbfile)
                except:
                    logging.warn("could not remove tmp file: '%s'", dbfile)

                if self.__background_file == None:
                    try:
                        os.remove(bgfile)
                    except:
                        logging.warn("could not remove tmp file: '%s'", bgfile)


    def make_sequence_file(self, seqs):
        """Creates a FASTA file from a list of(feature_id, sequence)
        pairs"""
        filename = None
        with tempfile.NamedTemporaryFile(prefix='memeseqs',
                                         delete=False) as outfile:
            filename = outfile.name
            st.write_sequences_to_fasta_file(outfile, seqs)
        return filename

    def dust(self, fasta_file_path):  # pylint: disable-msg=R0201
        """runs the dust command on the specified FASTA file and
        returns a list of sequences. It is assumed that dust has
        a very simple interface: FASTA in, output on stdout"""
        output = subprocess.check_output(['dust', fasta_file_path])
        return output

    # pylint: disable-msg=W0613,R0201
    def meme(self, infile_path, bgfile_path, num_motifs, pspfile_path=None):
        """Please implement me"""
        logging.error("MemeSuite.meme() - please implement me")

    def mast(self, meme_outfile_path, database_file_path,
             bgfile_path):  # pylint: disable-msg=R0201
        """Please implement me"""
        logging.error("MemeSuite.mast() - please implement me")


class MemeSuite430(MemeSuite):
    """Version 4.3.0 of MEME"""

    def meme(self, infile_path, bgfile_path, num_motifs,
             previous_motif_infos=None, pspfile_path=None):
        """runs the meme command on the specified input file, background file
        and positional priors file. Returns a tuple of
        (list of MemeMotifInfo objects, meme output)
        """
        command = ['meme', infile_path, '-bfile', bgfile_path,
                   '-time', '600', '-dna', '-revcomp',
                   '-maxsize', '9999999', '-nmotifs', str(num_motifs),
                   '-evt', '1e9', '-minw', '6', '-maxw', str(self.max_width()),
                   '-mod',  'zoops', '-nostatus', '-text']
        # if determine the seed sequence (-cons parameter) for this MEME run
        # uses the PSSM with the smallest score that has an e-value lower
        # than 0.1
        if previous_motif_infos != None:
            max_evalue = 0.1
            min_evalue = 10000000.0
            min_motif_info = None
            for motif_info in previous_motif_infos:
                if motif_info.evalue < min_evalue:
                    min_evalue = motif_info.evalue
                    min_motif_info = motif_info
            if min_motif_info != None and min_motif_info.evalue < max_evalue:
                cons = min_motif_info.consensus_string().upper()
                logging.info("seeding MEME with good motif %s", cons)
                command.extend(['-cons', cons])

        if pspfile_path:
            command.extend(['-psp', pspfile_path])

        #logging.info("running: %s", " ".join(command))
        output = subprocess.check_output(command)
        return (read_meme_output(output, num_motifs), output)

    def mast(self, meme_outfile_path, database_file_path,
             bgfile_path):
        """runs the mast command"""
        # note: originally run with -ev 99999, but MAST will crash with
        # memory errors
        command = ['mast', meme_outfile_path, '-d', database_file_path,
                   '-bfile', bgfile_path, '-nostatus', '-stdout', '-text',
                   '-brief', '-ev', '99999', '-mev', '99999', '-mt', '0.99',
                   '-seqp', '-remcorr']
        output = subprocess.check_output(command)
        return output


class MemeSuite481(MemeSuite):
    """Version 4.8.1 of MEME"""

    def meme(self, infile_path, bgfile_path, num_motifs, pspfile_path=None):
        """runs the meme command on the specified input file, background file
        and positional priors file. Returns a tuple of
        (list of MemeMotifInfo objects, meme output)
        """
        command = ['meme', infile_path, '-bfile', bgfile_path,
                   '-time', '600', '-dna', '-revcomp',
                   '-maxsize', '9999999', '-nmotifs', str(num_motifs),
                   '-evt', '1e9', '-minw', '6', '-maxw', str(self.max_width()),
                   '-mod',  'zoops', '-nostatus', '-text']

        if pspfile_path:
            command.append(['-psp', pspfile_path])

        #logging.info("running: %s", " ".join(command))
        output = subprocess.check_output(command)
        return (read_meme_output(output, num_motifs), output)

    def mast(self, meme_outfile_path, database_file_path,
             bgfile_path):
        """runs the mast command"""
        # note: originally run with -ev 99999, but MAST will crash with
        # memory errors
        command = ['mast', meme_outfile_path, database_file_path,
                   '-bfile', bgfile_path, '-nostatus', '-hit_list',
                   '-ev', '1500', '-mev', '99999', '-mt', '0.99',
                   '-seqp', '-remcorr']
        #logging.info("running: %s", " ".join(command))
        output = subprocess.check_output(command)
        return output


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


def read_meme_output(output_text, num_motifs):
    """Reads meme output file into a list of MotifInfo objects"""

    def extract_width(infoline):
        """extract the width value from the info line"""
        return int(__extract_regex('width =\s+\d+', infoline))

    def extract_num_sites(infoline):
        """extract the sites value from the info line"""
        return int(__extract_regex('sites =\s+\d+', infoline))

    def extract_llr(infoline):
        """extract the llr value from the info line"""
        return int(__extract_regex('llr =\s+\d+', infoline))

    def extract_evalue(infoline):
        """extract the e-value from the info line"""
        return float(__extract_regex('E-value =\s+\S+', infoline))

    def next_info_line(motif_number, lines):
        """finds the index of the next info line for the specified motif number
        1-based """
        return __next_regex_index('MOTIF\s+' + str(motif_number) + '.*',
                                  0, lines)

    def next_sites_index(start_index, lines):
        """returns the next sites index"""
        return __next_regex_index('[\t]Motif \d+ sites sorted by position ' +
                                  'p-value', start_index, lines)

    def read_sites(start_index, lines):
        """reads the sites"""
        sites_index = next_sites_index(start_index, lines)
        pattern = re.compile(
            "(\S+)\s+([+-])\s+(\d+)\s+(\S+)\s+\S+ (\S+) (\S+)?")
        current_index = sites_index + 4
        line = lines[current_index]
        sites = []
        while not line.startswith('----------------------'):
            match = pattern.match(line)
            sites.append((match.group(1), match.group(2), int(match.group(3)),
                          float(match.group(4)), match.group(5)))
            current_index += 1
            line = lines[current_index]
        return sites

    def read_pssm(start_index, lines):
        """reads the PSSM, in this case it's what is called the probability
        matrix in the meme output"""
        pssm_index = next_pssm_index(start_index, lines)
        current_index = pssm_index + 3
        line = lines[current_index]
        pattern = re.compile("\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")
        rows = []
        while not line.startswith('----------------------'):
            match = pattern.match(line)
            rows.append([float(match.group(1)), float(match.group(2)),
                         float(match.group(3)), float(match.group(4))])
            current_index += 1
            line = lines[current_index]
        return rows

    def next_pssm_index(start_index, lines):
        """determines the next PSSM start index"""
        return __next_regex_index('[\t]Motif \d+ position-specific ' +
                                  'probability matrix', start_index, lines)

    def read_motif_info(motif_number, lines):
        """Reads the MemeMotifInfo with the specified number from the input"""
        info_line_index = next_info_line(motif_number, lines)
        info_line = lines[info_line_index]
        return MemeMotifInfo(read_pssm(info_line_index + 1, lines),
                             motif_number,
                             extract_width(info_line),
                             extract_num_sites(info_line),
                             extract_llr(info_line),
                             extract_evalue(info_line),
                             read_sites(info_line_index + 1, lines))

    lines = output_text.split('\n')
    result = []
    for motif_number in xrange(1, num_motifs + 1):
        result.append(read_motif_info(motif_number, lines))
    return result


def read_mast_output(output_text, genes):
    """Reads out the p-values and e-values and the gene annotations
    from a mast output file"""
    def next_pe_value_line(start_index, lines):
        """Find the next combined p-value and e-value line"""
        return __next_regex_index('.*COMBINED P-VALUE.*',
                                  start_index, lines)

    def read_pe_values(lines):
        """read all combined p-values and e-values"""
        result = []
        current_index = next_pe_value_line(0, lines)
        while current_index != -1:
            gene = lines[current_index - 2].strip()
            line = lines[current_index]
            pvalue = float(__extract_regex('P-VALUE\s+=\s+(\S+)', line))
            evalue = float(__extract_regex('E-VALUE\s+=\s+(\S+)', line))
            result.append((gene, pvalue, evalue))
            current_index = next_pe_value_line(current_index + 1, lines)
        return result

    def read_seqalign_blocks(lines, start_index, seqlen):
        """Read the sequence alignment blocks starting at start_index
        a block has the format:
        1. motif number line (+/- = forward/reverse)
        2. pvalue line
        3. motif sequence line
        4. alignment/match line
        5. gene sequence line
        6. blank line (separator)
        -> Repeat this pattern until the whole database sequence printed

        While the mast output is easily human-readable, it
        is hard to parse programmatically.
        This method does it as follows:
        - read all motif numbers in sequence
        - read all p-values in sequencs
        - the motif number opening brackets are regarded as position markers

        for each block, we only need to keep track in which column the gene
        sequence starts and at which relative position we are
        """
        current_index = start_index
        is_last = False

        # global lines
        motifnum_line = ""
        pvalue_line = ""
        seq_line = ""
        while not is_last:
            is_last = is_last_block(lines, current_index, seqlen)
            # append to the motif number line, the p-value line, and seq line
            motifnum_line += lines[current_index].rstrip().ljust(80)[5:]
            pvalue_line += lines[current_index + 1].rstrip().ljust(80)[5:]
            seq_line += lines[current_index + 4].rstrip().ljust(80)[5:]
            current_index += 6

        motif_nums = read_motif_numbers(motifnum_line)
        positions = read_positions(motifnum_line, seq_line)
        pvalues = read_pvalues(pvalue_line, [(pos - 2) for pos in positions])
        return zip(pvalues, positions, motif_nums)

    def read_motifnum_line(line):
        """format and pad a motif number line"""
        return line.rstrip().ljust(80)[5:]

    def is_last_block(lines, index, seqlen):
        """determines whether the specified block is the last one for
        the current gene"""
        seqline = None
        try:
            seqline = lines[index + 4]
            seqstart_index = int(re.match('(\d+).*', seqline).group(1))
            seq_start = re.match('\d+\s+(\S+)', seqline).start(1)
            return ((len(seqline) - seq_start) + seqstart_index >= seqlen or
                    not re.match('(\d+).*', lines[index + 10]))
        except:
            if seqline != None:
                print "ERROR IN SEQLINE: [%s]" % seqline

    def read_motif_numbers(motifnum_line):
        """reads the motif numbers contained in a motif number line"""
        return [int(re.sub('\[|\]', '', motifnum))
                for motifnum in re.split(' +', motifnum_line)
                if len(motifnum.strip()) > 0]

    def read_pvalues(pvalue_line, indexes):
        """reads the p-values contained in a p-value line"""
        pvalues = []
        for index_num in xrange(len(indexes)):
            if index_num < len(indexes) - 1:
                pvalues.append(
                    float(pvalue_line[indexes[index_num]:
                                          indexes[index_num + 1]]))
            else:
                pvalues.append(float(
                        pvalue_line[indexes[index_num]:]))
        return pvalues

    def read_positions(motifnum_line, seqline):
        """we only need the motif number line and the sequence line
        to retrieve the position"""
        # offset +2 for compatibility with cMonkey R, don't really
        # know why we need this
        try:
            return [(m.start() + 2)
                    for m in re.finditer('\[', motifnum_line)]
        except:
            logging.error("ERROR in read_positions(), motifnum_line: '%s'",
                          str(motifnum_line))

    def read_annotations(lines, genes):
        """extract annotations, genes are given as refseq ids"""
        result = {}
        current_index = next_pe_value_line(0, lines)
        while current_index != -1:
            gene = lines[current_index - 2].strip()
            if gene in genes:
                info_line = lines[current_index]
                length = int(__extract_regex('LENGTH\s+=\s+(\d+)', info_line))
                has_seqalign_block = True
                diagram_match = re.match('^\s+DIAGRAM:\s+(\d+)$',
                                         lines[current_index + 1])
                if diagram_match != None:
                    diagram = int(diagram_match.group(1))
                    if diagram == length:
                        has_seqalign_block = False

                if has_seqalign_block:
                    # the diagram line can span several lines and the blank
                    # line after those can span several, so search for the
                    # first non-blank line after the block of blank lines
                    blank_index = current_index + 2
                    while len(lines[blank_index].strip()) > 0:
                        blank_index += 1
                    non_blank_index = blank_index + 1
                    while len(lines[non_blank_index].strip()) == 0:
                        non_blank_index += 1
                    result[gene] = read_seqalign_blocks(lines,
                                                        non_blank_index,
                                                        length)

            current_index = next_pe_value_line(current_index + 1, lines)
        return result

    # Make sure MAST returns a meaningful result
    if output_text.startswith("Error reading log-odds matrix file"):
        logging.warn("MAST returned the famous 'Error reading log-odds " +
                     "matrix file, provide empty result...'")
        return ([], {})
    else:
        lines = output_text.split('\n')
        pe_values = read_pe_values(lines)
        annotations = read_annotations(lines, genes)
        return (pe_values, annotations)


# extraction helpers
def __extract_regex(pattern, infoline):
    """generic info line field extraction based on regex"""
    try:
        match = re.search(pattern, infoline)
        return infoline[match.start():match.end()].split('=')[1].strip()
    except:
        logging.error("ERROR in __extract_regex(), pattern: '%s', infoline: '%s'",
                      str(pattern), str(infoline))


def __next_regex_index(pat, start_index, lines):
    """finds the line index of the first occurrence of the pattern"""
    line_index = start_index
    pattern = re.compile(pat)
    current_line = lines[line_index]
    while not pattern.match(current_line):
        line_index += 1
        if line_index >= len(lines):
            return -1
        current_line = lines[line_index]
    return line_index


def make_background_file(bgseqs, use_revcomp, bgorder=3):
    """create a meme background file and returns its name"""
    def make_seqs(seqs):
        """prepare the input sequences for feeding into meme.
        This means only taking the unique sequences and their reverse
        complement if desired"""
        meme_input_seqs = []
        for locseq in seqs.values():
            seq = locseq[1]
            util.add_if_unique(meme_input_seqs, seq)
            if use_revcomp:
                util.add_if_unique(meme_input_seqs, st.revcomp(seq))
        return meme_input_seqs

    filename = None
    bgmodel = st.markov_background(make_seqs(bgseqs), bgorder)
    with tempfile.NamedTemporaryFile(prefix='memebg',
                                     delete=False) as outfile:
        filename = outfile.name
        #logging.info("make background file '%s'", filename)
        outfile.write("# %s order Markov background model\n" %
                      util.order2string(len(bgmodel) - 1))
        for order_row in bgmodel:
            for seq, frequency in order_row.items():
                outfile.write('%s %10s\n' %
                              (seq, str(round(frequency, 8))))
    return filename


def global_background_file(organism, gene_aliases, seqtype, bgorder=3,
                           use_revcomp=True):
    """returns a background file that was computed on the set of all
    used sequences"""
    global_seqs = organism.sequences_for_genes_scan(gene_aliases,
                                                    seqtype=seqtype)
    logging.info("Computing global background file on seqtype '%s' " +
                 "(%d sequences)", seqtype, len(global_seqs))
    return make_background_file(global_seqs, use_revcomp, bgorder)

__all__ = ['read_meme_output']
