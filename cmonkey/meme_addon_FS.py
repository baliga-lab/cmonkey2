'''
Created on Apr 2, 2012

@author: frank
'''

import subprocess
import tempfile
import logging
import seqtools as st
import os
import util
import re
import numpy as np
import cElementTree as cET


MEME_CUTOFF = '1e9'     # lazy value to pass ALL

MAST_MIN_EVALUE = 1500 # original 1500
MAST_MEV = 99999




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
            remove_tempfiles=True, resultsDir = "meme_dummy"):
        """Create MemeSuite instance"""
        self._max_width = max_width
        self._use_revcomp = use_revcomp
        self._background_file = background_file
        self._remove_tempfiles = remove_tempfiles
        self._resultsDir = resultsDir

    def global_background_file(self):
        """returns the global background file used with this meme suite
        instance"""
        return self._background_file

    def max_width(self):
        """returns the max_width attribute"""
        return self._max_width

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
            if len(seq[1]) > self._max_width:
                seqs_for_dust[feature_id] = seq
        # only non-empty-input gets into dust, dust can not
        # handle empty input
        if len(seqs_for_dust) > 0:
            return process_with_dust(seqs_for_dust)
        else:
            return {}

    def __call__(self, input_seqs, all_seqs, iteration, cluster):
        return self.run_meme(input_seqs, all_seqs, iteration, cluster)

    def run_meme(self, input_seqs, all_seqs, iteration=0, cluster=0):
        """Runs the meme tool. input_seqs is a dictionary of
        (feature_id : (location, sequence)) that are to be provided as meme
        input, all_seqs is a dictionary that provides all sequences used
        in the cMonkey run, which will be used to compute background
        distribution"""
        def background_file():
            """decide whether to use global or specific background file"""
            if self._background_file != None:
                return self._background_file
            else:
                bgseqs = {feature_id: all_seqs[feature_id]
                          for feature_id in all_seqs
                          if feature_id not in input_seqs}
                logging.info("de novo creating background file")
                return make_background_file_FS(bgseqs, self._use_revcomp, "IamSupposedToBeATempFile", temp = True)

        #logging.info("run_meme() - iteration %d cluster %d # seqs = %d" % (iteration, cluster, len(input_seqs)))
        print "iteration %d cluster %d # seqs = %d" % (iteration, cluster, len(input_seqs)),
        bgfile = background_file()
        #logging.info("using background %s" %bgfile)
        seqfile = self.make_sequence_file(input_seqs.items())
        #logging.info("created sequence file in %s", seqfile)
        
        # Frank Schmitz
        # catch exceptions if meme fails:

        try:
            motif_infos, output = self.meme(seqfile, bgfile)
        except:
            logging.debug('meme failed iteration %d cluster %d # seqs = %d' % (iteration, cluster, len(input_seqs)))
            return MemeRunResult([], {}, [])
        #
        #
        #
        
        
        # run mast
        
        # generate dir for meme results
        output_dir_iter_base = "/".join([self._resultsDir, "Iteration%05d"%(iteration)])
        if not os.path.exists(output_dir_iter_base):
            os.mkdir(output_dir_iter_base)

        output_dir_meme_base = "/".join([output_dir_iter_base, "meme"])
        if not os.path.exists(output_dir_meme_base):
            os.mkdir(output_dir_meme_base)
        output_dir_mast_base = "/".join([output_dir_iter_base, "mast"])
        if not os.path.exists(output_dir_mast_base):
            os.mkdir(output_dir_mast_base)

        meme_outfile = "/".join([output_dir_meme_base, "cluster%04d_meme.txt"%(cluster)])
        with open(meme_outfile, "w") as outfile:
            outfile.write(output)
            outfile.close()
        #logging.info('wrote meme output to %s', meme_outfile)
        dbfile = self.make_sequence_file(
            [(feature_id, locseq[1])
             for feature_id, locseq in all_seqs.items()])
        #logging.info('created mast database in %s', dbfile)



        output_dir_mast = "/".join([output_dir_mast_base, "cluster%04d"%(cluster)])
        if not os.path.exists(output_dir_mast):
            os.mkdir(output_dir_mast)
        try:
            mastfile, output = self.mast(meme_outfile, dbfile, bgfile, output_dir_mast)
            pe_values, annotations = read_mast_output(mastfile, output,
                                                      input_seqs.keys())

            return MemeRunResult(pe_values, annotations, motif_infos)
        except:
            logging.info('read mast output failed, no returns here')
            return MemeRunResult([], {}, [])
        finally:
            if self._remove_tempfiles:
                logging.info("DELETING ALL TMP FILES...")
                
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
    
                #if self._background_file == None:
                #    try:
                #        os.remove(bgfile)
                #    except:
                #        logging.warn("could not remove tmp file: '%s'", bgfile)


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
    def meme(self, infile_path, bgfile_path, num_motifs=2,
             pspfile_path=None):
        """Please implement me"""
        logging.error("MemeSuite.meme() - please implement me")

    def mast(self, meme_outfile_path, database_file_path,
             bgfile_path):  # pylint: disable-msg=R0201
        """Please implement me"""
        logging.error("MemeSuite.mast() - please implement me")


class MemeSuite481(MemeSuite):
    """Version 4.8.1 of MEME"""

    def __init__(self, max_width=24, use_revcomp=True, background_file=None,
            remove_tempfiles=True, resultsDir = "meme_dummy"):
        self._num_motifs = 2
        
        MemeSuite.__init__(self, max_width = max_width,
                           use_revcomp = use_revcomp,
                           background_file = background_file,
                           remove_tempfiles = remove_tempfiles,
                           resultsDir = resultsDir)

    def meme(self, infile_path, bgfile_path, pspfile_path=None):
        """runs the meme command on the specified input file, background file
        and positional priors file. Returns a tuple of
        (list of MemeMotifInfo objects, meme output)
        """
        command = ['meme', infile_path, '-bfile', bgfile_path,
                   '-time', '600', '-dna', '-revcomp',
                   '-maxsize', '9999999', '-nmotifs', str(self._num_motifs),
                   '-evt', '1e9', '-minw', '6', '-maxw', str(self.max_width()),
                   '-mod',  'zoops', '-nostatus', '-text']
        if pspfile_path:
            command.append(['-psp', pspfile_path])

        #logging.info("running: %s", " ".join(command))
        output = subprocess.check_output(command)
        return (read_meme_output(output, self._num_motifs), output)

    def mast(self, meme_outfile_path, database_file_path,
             bgfile_path, output_dir):
        """runs the mast command"""
        # note: originally run with -ev 99999, but MAST will crash with
        # memory errors
        command = ['mast', meme_outfile_path, database_file_path,
                   '-bfile', bgfile_path, '-nostatus', '-oc', output_dir,
                   '-ev', '1500', '-mev', '99999', '-mt', '0.99',
                   '-seqp', '-remcorr', '-notext', '-nohtml']
        #logging.info("running: %s", " ".join(command))
        output = subprocess.check_output(command)
        filename = "/".join([output_dir, "mast.xml"])
        return (filename, output)


class MemeMotifInfo:
    """Only a motif's info line, the
    probability matrix and the site information is relevant"""
    # pylint: disable-msg=R0913
    def __init__(self, pssm, motif_num, width, num_sites, llr, evalue, sites):
        """Creates a MemeMotifInfo instance"""
        self.__pssm = pssm
        self.__motif_num = motif_num
        self.__width = width
        self.__num_sites = num_sites
        self.__llr = llr
        self.__evalue = evalue
        self.__sites = sites


    def motif_num(self):
        """returns the motif number"""
        return self.__motif_num

    def pssm(self):
        """return the PSSM rows"""
        return self.__pssm

    def width(self):
        """Returns the width"""
        return self.__width

    def num_sites(self):
        """returns the number of sites"""
        return self.__num_sites

    def llr(self):
        """returns the log-likelihood ratio"""
        return self.__llr

    def evalue(self):
        """returns the e value"""
        return self.__evalue

    def sites(self):
        """returns the sites"""
        return self.__sites

    def consensus_string(self, cutoff1=0.7, cutoff2=0.4):
        """returns the consensus string from the pssm table
        remember: letter order is ACGT"""
        alphabet = 'ACGT'
        result = ""
        for row in xrange(len(self.__pssm)):
            rowvals = self.__pssm[row]
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
        return ("Motif width: %d sites: %d llr: %d e-value: %f" %
         (self.width(), self.num_sites(), self.llr(),
          self.evalue()))


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



class MastResult():
    def __init__(self, genes):
        self.__command_line = None
        self.__max_correlation = None
        self.__remove_correlated = None
        self.__strand_handling = None
        self.__translate_dna = None
        self.__max_seq_evalue = None
        self.__adj_hit_pvalue = None
        self.__max_hit_pvalue = None
        self.__max_weak_pvalue = None
        self.__host = None
        self.__when = None
        
        self.__valueDict = {}
        self.__annotDict = {}

        self.__pe_values = []
        
        self.__genes = dict(zip(genes, genes))  # funny, but makes searching MUCH faster
        
    def model(self, elem):
        """ parse the model part of the xml file
        den teil kann man sich wohl auch sparen
        
        """
        for subelem in elem:
            if subelem == "command_line":
                self.__command_line = subelem.text
            elif subelem == "max_correlation":
                self.__max_correlation = float(subelem.text)
            elif subelem == "max_seq_evalue":
                self.__max_seq_evalue = float(subelem.text)
            elif subelem == "max_hit_pvalue":
                self.__max_hit_pvalue = float(subelem.text)
            elif subelem == "max_weak_pvalue":
                self.__max_weak_pvalue = float(subelem.text)
            elif subelem == "host":
                self.__host = subelem.text
            elif subelem == "when":
                self.__when = subelem.text
                
            elif subelem == "remove_correlated":
                for name, value in subelem.items():
                    if name == 'value': self.__remove_correlated = value
            elif subelem == "strand_handling":
                for name, value in subelem.items():
                    if name == 'value': self.__strand_handling = value
            elif subelem == "translate_dna":
                for name, value in subelem.items():
                    if name == 'value': self.__translate_dna = value
            elif subelem == "adj_hit_pvalue":
                for name, value in subelem.items():
                    if name == 'value': self.__adj_hit_pvalue = value

            
    def alphabet(self, elem):
        """
        not implemented yet, if needed at all??
        """
        pass
    def motifs(self, elem):
        """
        not implemented yet, if needed at all??
        """
        pass
    def sequences(self, elem):
        def workoutseq(seqelem):
            hitDict = {}
            annotL = []
            
            for tup in seqelem.items():
                if tup[0] == 'name': seqName = tup[1]
            
            # find the score line
            subby = seqelem.find('score')
            for name, value in subby.items():
                if name == 'combined_pvalue': cp_value = float(value)
                elif name == 'evalue': evalue = float(value)
            
            # here, the info we need is all, we can return and save time
            if seqName not in self.__genes: return (seqName, cp_value, evalue, annotL)
            # if not, find all the hits with scores and build the annotation List 
            for subelem in seqelem:
                if subelem.tag != 'seg': continue
                for subsubelem in subelem:
                    if subsubelem.tag != 'hit': continue
                    hitDict = dict(subsubelem.items())
                    if hitDict['strand'] == 'reverse':
                        motif = (-1) * int (hitDict['motif'].replace('motif_', ''))
                    elif hitDict['strand'] == 'forward':
                        motif = int (hitDict['motif'].replace('motif_', ''))
                    annotL.append((float(hitDict['pvalue']),
                                  int(hitDict['pos']),
                                  motif ))
            return (seqName, cp_value, evalue, annotL)
                

        
        for Sequence in elem:
            if Sequence.tag == "database":
                for name, value in Sequence.items():
                    if name == 'name': self.__dbname = value
                    elif name == 'seq_count': self.__seq_count = int(value) 
            elif Sequence.tag == "sequence":
                SeqName, cp_value, evalue, hitList = workoutseq(Sequence)
                self.__valueDict[SeqName] = (cp_value, evalue)
                if SeqName in self.__genes:
                    self.__annotDict[SeqName] = hitList
    
    
    def runtime(self, elem):
        """
        not implemented yet, if needed at all??
        """
        pass
    

    def getSeqValueD(self):
        return self.__valueDict

    def get_pe_values(self):
        return self.__pe_values

    def getAnnotD(self):
        return self.__annotDict

    def combine_p_values(self):
        """ the frankish way of combining multiple sequence scores for many genes
        """
        
        geneL = []
        pe_valL = []
        gpvD = {}
        gEsD = {}
        for geneID, scoretup in self.__valueDict.iteritems():
            pval = scoretup[0]
            Escore = scoretup[1]
            gene = geneID.split("_")[0]
            try:
                temppvL = gpvD[gene]
                tempEsL = gEsD[gene]
            except KeyError:
                temppvL, tempEsL = [], []
                geneL.append(gene)
            temppvL.append(pval)
            tempEsL.append(Escore)
            gpvD[gene] = temppvL
            gEsD[gene] = tempEsL
        for gene in geneL:
            pvmin = np.array(gpvD[gene]).min()
            Esmin = np.array(gEsD[gene]).max()
            pe_valL.append((gene, pvmin, Esmin))
        self.__pe_values = pe_valL

        




def read_mast_output(mastfile, output, genes):
    MAST = MastResult(genes)
    for event, elem in cET.iterparse(mastfile):
        tag = elem.tag
        if tag == "model":  
            MAST.model(elem)
        elif tag == "sequences":  
            MAST.sequences(elem)
    annotations = MAST.getAnnotD()
    MAST.combine_p_values()
    comb_pe_values = MAST.get_pe_values() 
    return (comb_pe_values, annotations)




# extraction helpers
def __extract_regex(pattern, infoline):
    """generic info line field extraction based on regex"""
    match = re.search(pattern, infoline)
    return infoline[match.start():match.end()].split('=')[1].strip()


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




# Frank Schmitz:
# a new variant of the make background file
# puts all seqs into a gigantic string, divided by spaces
# this one can be searched through regex without iterations
# regexs are a child of the devil...hell they're fast
# but - in the case of repeat masked files, it is till terribly slow
# consider another alternative to 'fill in' the masked

def make_background_file_FS(bgseqs, use_revcomp, filename, alphabet, alphabet_replacement, temp = True):
    """create a meme background file and returns its name"""
    logging.info("\x1b[31mmeme:\t\x1b[0mgenerating Seq background model")
    
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
        bigstring = "  ".join(meme_input_seqs)
        return bigstring

    bgmodel = st.markov_background_FS(make_seqs(bgseqs), 3, alphabet, alphabet_replacement)
    logging.info("\x1b[31mmeme:\t\x1b[0mMarkov model built")

    if temp == False:
        outfile = open(filename, "w")
    else:
        filename = None
        outfile = tempfile.NamedTemporaryFile(prefix='memebg',
                                         delete=False)
        filename = outfile.name

    
    logging.info("\x1b[31mmeme:\t\x1b[0mmake background file '%s'", filename)
    outfile.write("# %s order Markov background model\n" %
                  util.order2string(len(bgmodel) - 1))
    for order_row in bgmodel:
        for seq, frequency in order_row.items():
            outfile.write('%s %10s\n' %
                          (seq, str(round(frequency, 8))))
    if temp == False: outfile.close()
    logging.info("\x1b[31mmeme:\t\x1b[0mdone generating background Seqs")
    return filename





def global_background_file_FS(organism, gene_aliases, seqtype, alphabet, alphabet_replacement, use_revcomp=True, resultsDir = ""):
    """returns a background file that was computed on the set of all
    used sequences"""
    global_seqs = organism.sequences_for_genes_scan(gene_aliases,
                                                    seqtype=seqtype)
    logging.info("\x1b[31mmeme:\t\x1b[0mComputing global background file on seqtype '%s' " +
                 "(%d sequences)", seqtype, len(global_seqs))
    if not os.path.exists(resultsDir):
        os.mkdir(resultsDir)
    filename = "/".join([resultsDir, "memebg.txt"])

    return make_background_file_FS(global_seqs, use_revcomp, filename, alphabet, alphabet_replacement, temp = False)

