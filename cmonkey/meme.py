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
import shutil
import re
import collections
import xml.etree.ElementTree as ET


MemeRunResult = collections.namedtuple('MemeRunResult',
                                       ['pe_values', 'annotations', 'motif_infos'])


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
    def __init__(self, config_params, background_file=None, bgmodel=None,
                 remove_tempfiles=True):
        """Create MemeSuite instance"""
        self.max_width = int(config_params['MEME']['max_width'])
        self.background_order = int(config_params['MEME']['background_order'])
        self.__use_revcomp = config_params['MEME']['use_revcomp'] == 'True'
        self.__background_file = background_file
        self.bgmodel = bgmodel
        self.__remove_tempfiles = remove_tempfiles
        self.arg_mod = config_params['MEME']['arg_mod']

    def global_background_file(self):
        """returns the global background file used with this meme suite
        instance"""
        return self.__background_file

    def remove_low_complexity(self, seqs):
        """send sequences through dust filter, send only those
        to dust that are larger than max_width"""
        def process_with_dust(seqs):
            """data conversion from and to dust tool"""
            dust_tmp_file = None
            with tempfile.NamedTemporaryFile(prefix='dust',
                                             delete=False) as dust_input:
                for feature_id, seq in seqs.iteritems():
                    dust_input.write(">%s\n" % feature_id)
                    if isinstance(seq, str):
                        dust_input.write("%s\n" % seq)
                    else:
                        dust_input.write("%s\n" % seq)
                dust_tmp_file = dust_input.name
                #logging.info("DUST input written to: %s", dust_input.name)
            seqpairs = st.read_sequences_from_fasta_string(
                self.dust(dust_tmp_file))
            os.remove(dust_tmp_file)
            return {feature_id: seq for feature_id, seq in seqpairs}

        seqs_for_dust = {}
        for feature_id, seq in seqs.iteritems():
            if isinstance(seq, str):
                if len(seq) > self.max_width:
                    seqs_for_dust[feature_id] = seq
            else:
                if len(seq[1]) > self.max_width:
                    seqs_for_dust[feature_id] = seq[1]
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
        feature_ids = set(params.feature_ids)  # optimization: reduce lookup time
        input_seqs = params.seqs
        all_seqs = params.used_seqs

        def background_file():
            """decide whether to use global or specific background file"""
            if self.__background_file is not None:
                #logging.info("using global background: '%s'", self.__background_file)
                return self.__background_file
            else:
                bgseqs = {feature_id: all_seqs[feature_id]
                          for feature_id in all_seqs
                          if feature_id not in feature_ids}
                # note: the result of make_background_file is a tuple !!
                return make_background_file(bgseqs, self.__use_revcomp,
                                            self.background_order)[0]

        #logging.info("run_meme() - # seqs = %d", len(input_seqs))
        bgfile = background_file()
        #logging.info("created background file in %s", bgfile)
        seqfile = self.make_sequence_file(
            [(feature_id, input_seqs[feature_id])
             for feature_id in params.feature_ids if feature_id in input_seqs])
        #logging.info("created sequence file in %s", seqfile)
        motif_infos, output = self.meme(seqfile, bgfile, params.num_motifs,
                                        previous_motif_infos=params.previous_motif_infos)

        # run mast
        meme_outfile = None
        is_last_iteration = params.iteration > params.num_iterations
        if 'keep_memeout' in params.debug or is_last_iteration:
            meme_outfile = os.path.join(params.outdir,
                                        'meme-out-%04d-%04d' % (params.iteration, params.cluster))
            with open(meme_outfile, 'w') as outfile:
                outfile.write(output)
        else:
            with tempfile.NamedTemporaryFile(prefix='meme.out.',
                                             delete=False) as outfile:
                meme_outfile = outfile.name
                outfile.write(output)

        #logging.info('wrote meme output to %s', meme_outfile)
        dbfile = self.make_sequence_file(
            [(feature_id, locseq[1])
             for feature_id, locseq in all_seqs.iteritems()])
        #logging.info('created mast database in %s', dbfile)
        try:
            mast_output = self.mast(meme_outfile, dbfile, bgfile)
            if 'keep_mastout' in params.debug:
                with open('%s.mast' % meme_outfile, 'w') as outfile:
                    outfile.write(mast_output)
            pe_values, annotations = self.read_mast_output(mast_output,
                                                           input_seqs.keys())
            return MemeRunResult(pe_values, annotations, motif_infos)
        except subprocess.CalledProcessError, e:
            if e.output.startswith('No input motifs pass the E-value'):
                logging.warn("no input motifs pass the e-value, ignoring result")
                return MemeRunResult([], [], [])
            else:
                print "Unknown error in MAST:\n ", e.__dict__
                logging.error("MAST error: %s", e.output)
                raise
        finally:
            if self.__remove_tempfiles:
                try:
                    os.remove(seqfile)
                except:
                    logging.warn("could not remove tmp file: '%s'", seqfile)
                try:
                    if 'keep_memeout' not in params.debug and not is_last_iteration:
                        os.remove(meme_outfile)
                except:
                    logging.warn("could not remove tmp file: '%s'", meme_outfile)
                try:
                    os.remove(dbfile)
                except:
                    logging.warn("could not remove tmp file: '%s'", dbfile)

                if self.__background_file is None:
                    try:
                        os.remove(bgfile)
                    except:
                        logging.warn("could not remove tmp file: '%s'", bgfile)

    def read_mast_output(self, mast_output, genes):
        """Please implement me"""
        logging.error("MemeSuite.read_mast_output() - please implement me")

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
    def meme(self, infile_path, bgfile_path, num_motifs,
             previous_motif_infos=None, pspfile_path=None):
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
                   '-evt', '1e9', '-minw', '6', '-maxw', str(self.max_width),
                   '-mod',  self.arg_mod, '-nostatus', '-text']
        # if determine the seed sequence (-cons parameter) for this MEME run
        # uses the PSSM with the smallest score that has an e-value lower
        # than 0.1
        if previous_motif_infos is not None:
            max_evalue = 0.1
            min_evalue = 10000000.0
            min_motif_info = None
            for motif_info in previous_motif_infos:
                if motif_info.evalue < min_evalue:
                    min_evalue = motif_info.evalue
                    min_motif_info = motif_info
            if min_motif_info is not None and min_motif_info.evalue < max_evalue:
                cons = min_motif_info.consensus_string().upper()
                logging.debug("seeding MEME with good motif %s", cons)
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
                   '-brief', '-ev', '999999', '-mev', '9999999', '-mt', '0.99',
                   '-seqp', '-remcorr']
        #logging.info("running: %s", " ".join(command))
        output = subprocess.check_output(command, stderr=subprocess.STDOUT)
        return output

    def read_mast_output(self, mast_output, genes):
        """old-style MAST output"""
        return read_mast_output_oldstyle(mast_output, genes)


class MemeSuite481(MemeSuite):
    """Version 4.8.1 of MEME"""

    def meme(self, infile_path, bgfile_path, num_motifs,
             previous_motif_infos=None,
             pspfile_path=None):
        """runs the meme command on the specified input file, background file
        and positional priors file. Returns a tuple of
        (list of MemeMotifInfo objects, meme output)
        """
        command = ['meme', infile_path, '-bfile', bgfile_path,
                   '-time', '600', '-dna', '-revcomp',
                   '-maxsize', '9999999', '-nmotifs', str(num_motifs),
                   '-evt', '1e9', '-minw', '6', '-maxw', str(self.max_width),
                   '-mod',  self.arg_mod, '-nostatus', '-text']

        ### NOTE: There is a bug in current MEME 4.9.0, that can cause the
        ### -cons option to crash
        ### ----- We leave it out here for now
        """
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
        """
        if pspfile_path:
            command.extend(['-psp', pspfile_path])

        #logging.info("running: %s", " ".join(command))
        try:
            output = subprocess.check_output(command)
            return (read_meme_output(output, num_motifs), output)
        except:
            print command
            raise

    def mast(self, meme_outfile_path, database_file_path,
             bgfile_path):
        """runs the mast command. Version 4.81 and above behave differently
        than 4.30: The output will be generated in an output directory
        So, here we'll generate a temporary directory
        """
        # note: originally run with -ev 99999, but MAST will crash with
        # memory errors
        dirname = tempfile.mkdtemp(prefix="mastout")
        try:
            command = ['mast', meme_outfile_path, database_file_path,
                       '-bfile', bgfile_path, '-nostatus',
                       '-ev', '1500', '-mev', '99999', '-mt', '0.99', '-nohtml',
                       '-notext', '-seqp', '-remcorr', '-oc', dirname]
            logging.debug("running: %s", " ".join(command))
            output = subprocess.check_output(command, stderr=subprocess.STDOUT)
            with open(os.path.join(dirname, "mast.xml")) as infile:
                result = infile.read()
            return result
        except subprocess.CalledProcessError, e:
            logging.warn("there is an exception thrown in MAST: %s",
                         e.output)
            return None  # return nothing if there was an error
        finally:
            logging.debug("removing %s...", dirname)
            shutil.rmtree(dirname)
            ##print "done."

    def read_mast_output(self, mast_output, genes):
        """XML MAST output"""
        return read_mast_output_xml(mast_output, genes)


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
            if match is None:
                logging.error("ERROR in read_pssm(), line(#%d) is: '%s'", current_index, line)
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


def read_mast_output_xml(output_text, genes):
    """Reads p/e values and gene annotations from a MAST output file
    in XML format.
    Inputs: - output_text: a string in MAST XML output format
    ------- - genes: a list of genes that were used as input to
              the previous MEME run
    Returns: a pair (pevalues, annotations)
    -------- - pevalues is [(gene, pval, eval)]
             - annotations is a dictionary gene -> [(pval, pos, motifnum)]"""
    pevalues = []
    annotations = {}
    if output_text is None:  # there was an error in mast, ignore its output
        return pevalues, annotations

    root = ET.fromstring(output_text)
    for sequence in root.iter('sequence'):
        score = sequence.find('score')
        seqname = sequence.get('name')
        if not seqname in annotations:
            annotations[seqname] = []
        pevalues.append((seqname,
                         float(score.get('combined_pvalue')),
                         float(score.get('evalue'))))
        if seqname in genes:
            for hit in sequence.iter('hit'):
                strand = hit.get('strand')
                motifnum = int(hit.get('motif').replace('motif_', ''))
                if strand == 'reverse':
                    motifnum = -motifnum
                annot = (float(hit.get('pvalue')),
                         int(hit.get('pos')) + 2,  # like R cmonkey
                         motifnum)
                annotations[seqname].append(annot)
    return pevalues, annotations


def read_mast_output_oldstyle(output_text, genes):
    """Reads out the p-values and e-values and the gene annotations
    from a mast output file. This format is generated by
    MAST 4.30 and is only here to support the legacy format.
    Use the XML version instead, it is more reliable.
    """
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
            if seqline is not None:
                print "ERROR IN SEQLINE: [%s]" % seqline

    def read_motif_numbers(motifnum_line):
        """reads the motif numbers contained in a motif number line"""
        return [int(re.sub('\[|\]', '', motifnum))
                for motifnum in re.split(' +', motifnum_line)
                if len(motifnum.strip()) > 0]

    def read_pvalues(pvalue_line, indexes):
        """reads the p-values contained in a p-value line"""
        def make_float(s):
          """unfortunately, MEME result lines can have weird float formats"""
          return float(s.replace(' ', ''))
        pvalues = []
        for index_num in xrange(len(indexes)):
            if index_num < len(indexes) - 1:
                pvalues.append(
                    make_float(pvalue_line[indexes[index_num]:
                                      indexes[index_num + 1]]))
            else:
                pvalues.append(make_float(pvalue_line[indexes[index_num]:]))
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
                if diagram_match is not None:
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


def make_background_file(bgseqs, use_revcomp, bgorder):
    """create a meme background file and returns its name and the model itself as
    a tuple"""
    def make_seqs(seqs):
        """prepare the input sequences for feeding into meme.
        This means only taking the unique sequences and their reverse
        complement if desired"""
        meme_input_seqs = []
        for locseq in seqs.values():
            seq = locseq[1]
            if seq not in meme_input_seqs:
                meme_input_seqs.append(seq)
            if use_revcomp:
                revseq = st.revcomp(seq)
                if revseq not in meme_input_seqs:
                    meme_input_seqs.append(revseq)
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
            for seq, frequency in order_row.iteritems():
                outfile.write('%s %10s\n' %
                              (seq, str(round(frequency, 8))))
    return (filename, bgmodel)


def global_background_file(organism, gene_aliases, seqtype, bgorder=3,
                           use_revcomp=True):
    """returns a background file that was computed on the set of all
    used sequences"""
    global_seqs = organism.sequences_for_genes_scan(gene_aliases,
                                                    seqtype=seqtype)
    logging.debug("Computing global background file on seqtype '%s' " +
                  "(%d sequences)", seqtype, len(global_seqs))
    return make_background_file(global_seqs, use_revcomp, bgorder)


USER_TEST_FASTA_PATH = 'config/fasta_test.fa'
SYSTEM_TEST_FASTA_PATH = '/etc/cmonkey2/fasta_test.fa'


def check_meme_version():
    logging.info('checking MEME...')
    if os.path.exists(USER_TEST_FASTA_PATH):
        test_fasta = USER_TEST_FASTA_PATH
    elif os.path.exists(SYSTEM_TEST_FASTA_PATH):
        test_fasta = SYSTEM_TEST_FASTA_PATH
    else:
        raise Exception('fasta_test.fa not found !')

    try:
        command = ['meme', '-nostatus', '-text', test_fasta]
        output = subprocess.check_output(command).split('\n')
        for line in output:
            if line.startswith('MEME version'):
                return line.split(' ')[2]
    except OSError:
        logging.error("MEME does not exist")
        return None

######################################################################
### Export
###############################3

MEME_FILE_HEADER = """MEME version 3.0

ALPHABET= ACGT

strands: + -

Background letter frequencies (from dataset with add-one prior applied):
A %.3f C %.3f G %.3f T %.3f
"""

def write_pssm(outfile, cursor, motif_info_id, evalue, num_sites):
    """writes a single PSSM to the given file"""
    outfile.write('\nMOTIF %d\n' % motif_info_id)
    outfile.write('BL   MOTIF %s width=0 seqs=0\n' % motif_info_id)

    cursor.execute('select a,c,g,t from motif_pssm_rows where motif_info_id=? order by row',
                   [motif_info_id])
    pssm_rows = [(a, c, g, t) for a, c, g, t in cursor.fetchall()]
    outfile.write('letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e\n' % (len(pssm_rows), num_sites, evalue))
    for a, c, g, t in pssm_rows:
        outfile.write('%5.3f %5.3f %5.3f %5.3f\n' % (a, c, g, t))

def write_motifs2meme(conn, filepath):
    """Write the motifs to a MEME file and returns True if successful
    Currently, this only works if there is global background data in the database
    """
    cursor = conn.cursor()
    cursor2 = conn.cursor()

    cursor.execute("select subsequence,pvalue from global_background where subsequence in ('A','C','G','T')")
    freqs = {base: pvalue for base, pvalue in cursor.fetchall()}
    if len(freqs) >= 4 and 'A' in freqs and 'C' in freqs and 'G' in freqs and 'T' in freqs:
        logging.debug('retrieving letter frequency from global background distribution')
        with open(filepath, 'w') as outfile:
            outfile.write(MEME_FILE_HEADER % (freqs['A'], freqs['C'], freqs['G'], freqs['T']))
            cursor.execute('select max(iteration) from motif_infos')
            iteration = cursor.fetchone()[0]
            cursor.execute('select mi.rowid,evalue,count(mms.rowid) from motif_infos mi join meme_motif_sites mms on mi.rowid=mms.motif_info_id where iteration=? group by mi.rowid',
                           [iteration])
            num_pssms_written = 0
            for motif_info_id, evalue, num_sites in cursor.fetchall():
                write_pssm(outfile, cursor2, motif_info_id, evalue, num_sites)

            # no pssms were written, this can happen when we did not
            # run MEME, but weeder, so we retrieve the number of sites from MAST instead
            if num_pssms_written == 0:
                cursor.execute("select mi.rowid,evalue,count(ann.rowid) from motif_infos mi join motif_annotations ann on mi.rowid=ann.motif_info_id where mi.iteration=? group by mi.rowid",
                               [iteration])
                for motif_info_id, evalue, num_sites in cursor.fetchall():
                    write_pssm(outfile, cursor2, motif_info_id, evalue, num_sites)
        return True
    else:
        logging.warn('no global background distribution found')
        return False


EVALUE_CUTOFF = 100
RESID_CUTOFF  = 0.8
DIST_METHOD   = "ed"
Q_THRESHOLD   = 0.5
MIN_OVERLAP   = 4
Q_PSEUDO      = 0
T_PSEUDO      = 0


def run_tomtom(conn, targetdir, version, q_thresh=Q_THRESHOLD, dist_method=DIST_METHOD,
               min_overlap=MIN_OVERLAP, q_pseudo=Q_PSEUDO, t_pseudo=T_PSEUDO):
    """a wrapper around the tomtom script"""
    targetfile = os.path.join(targetdir, 'post.tomtom.meme')
    queryfile = targetfile
    if write_motifs2meme(conn, targetfile):
        command = ['tomtom',
                   '-verbosity', '1',
                   '-q-thresh', '%.3f' % q_thresh,
                   '-dist', dist_method,
                   '-min-overlap', '%d' % min_overlap,
                   '-text',
                   '-query-pseudo', '%.3f' % q_pseudo,
                   '-target-pseudo', '%.3f' % t_pseudo]
        logging.debug(" ".join(command))

        # Tomtom versions > 4.8.x drop the target and query switches
        if version == '4.3.0':
            command.extend(['-target', targetfile, '-query', queryfile])
        else:
            command.extend([targetfile, queryfile])

        try:
            output = subprocess.check_output(command)
            lines = output.split('\n')[1:]
            for line in lines:
                if len(line.strip()) > 0:
                    row = line.strip().split('\t')            
                    motif1 = int(row[0])
                    motif2 = int(row[1])
                    if motif1 != motif2:
                        pvalue = float(row[3])
                        conn.execute('insert into tomtom_results (motif_info_id1,motif_info_id2,pvalue) values (?,?,?)',
                                     [motif1, motif2, pvalue])
        except:
            raise

__all__ = ['read_meme_output']
