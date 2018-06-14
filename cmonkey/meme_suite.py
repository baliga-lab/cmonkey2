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
import os
import shutil
import re
import collections
import xml.etree.ElementTree as ET
from pkg_resources import Requirement, resource_filename, DistributionNotFound

import cmonkey.seqtools as st
import cmonkey.util as util
import cmonkey.database as cm2db

# For now until we integrate these better
import cmonkey.meme.meme as meme_formats
import cmonkey.meme.mast as mast_formats

from sqlalchemy import func

try:
    xrange
except NameError:
    xrange = range


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
            with tempfile.NamedTemporaryFile(mode='w+',
                                             prefix='dust',
                                             delete=False) as dust_input:
                for feature_id, seq in seqs.items():
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
        for feature_id, seq in seqs.items():
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

        try:
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
            mast_failed = False
            is_last_iteration = params.iteration > params.num_iterations
            if 'keep_memeout' in params.debug or is_last_iteration:
                meme_outfile = os.path.join(params.outdir,
                                            'meme-out-%04d-%04d' % (params.iteration, params.cluster))
                with open(meme_outfile, 'w') as outfile:
                    outfile.write(output)
            else:
                with tempfile.NamedTemporaryFile(mode='w+', prefix='meme.out.',
                                                 delete=False) as outfile:
                    meme_outfile = outfile.name
                    outfile.write(output)

            #logging.info('wrote meme output to %s', meme_outfile)
            dbfile = self.make_sequence_file(
                [(feature_id, locseq[1])
                 for feature_id, locseq in all_seqs.items()])
            #logging.info('created mast database in %s', dbfile)
        except subprocess.CalledProcessError as e:
            logging.error("MEME output: %s", e.output)
            return MemeRunResult([], [], [])

        try:
            mast_output = self.mast(meme_outfile, dbfile, bgfile)
            # There is a bug in MAST, catch that here to report to MEME team
            # when it is fixed, we could remove it
            if mast_output is None:
                mast_failed = True

            if 'keep_mastout' in params.debug:
                with open('%s.mast' % meme_outfile, 'w') as outfile:
                    outfile.write(mast_output)
            pe_values, annotations = self.read_mast_output(mast_output,
                                                           input_seqs.keys())
            return MemeRunResult(pe_values, annotations, motif_infos)
        except subprocess.CalledProcessError as e:
            if e.output.startswith('No input motifs pass the E-value'):
                logging.warn("no input motifs pass the e-value, ignoring result")
                return MemeRunResult([], [], [])
            else:
                print("Unknown error in MAST:\n ", e.__dict__)
                logging.error("MAST error: %s", e.output)
                return MemeRunResult([], [], [])
        finally:
            if mast_failed:
                # This is a workaround to keep the meme output file in case
                # the strange 1000 != 1000 text parser problem
                # occurs in MAST. Remove when it's fixed in MEME
                shutil.copyfile(seqfile, '/tmp/masterror-seqs')
                shutil.copyfile(meme_outfile, '/tmp/masterror-memeout')
                shutil.copyfile(dbfile, '/tmp/masterror-dbfile')
                shutil.copyfile(bgfile, '/tmp/masterror-bgfile')
            elif self.__remove_tempfiles:
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
        with tempfile.NamedTemporaryFile(mode='w+', prefix='memeseqs',
                                         delete=False) as outfile:
            filename = outfile.name
            st.write_sequences_to_fasta_file(outfile, seqs)
        return filename

    def dust(self, fasta_file_path):  # pylint: disable-msg=R0201
        """runs the dust command on the specified FASTA file and
        returns a list of sequences. It is assumed that dust has
        a very simple interface: FASTA in, output on stdout"""
        output = subprocess.check_output(['dust', fasta_file_path])
        return output.decode('utf-8')

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
        output = subprocess.check_output(command).decode('utf-8')
        return (meme_formats.from_test(output, num_motifs), output)

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
        return output.decode('utf-8')

    def read_mast_output(self, mast_output, genes):
        """old-style MAST output"""
        return mast_formats.from_430_text(mast_output, genes)


class MemeSuite481(MemeSuite):
    """Supports versions 4.8.1 and greater of MEME"""

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
        if previous_motif_infos is not None:
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
            output = subprocess.check_output(command).decode('utf-8')
            return (meme_formats.from_text(output, num_motifs), output)
        except:
            logging.error("MEME execution error, command: %s", str(command))
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
        except subprocess.CalledProcessError as e:
            logging.warn("there is an exception thrown in MAST: %s, (meme file: '%s', dbfile: '%s', bgfile: '%s')",
                         e.output, meme_outfile_path, database_file_path, bgfile_path)
            return None  # return nothing if there was an error
        finally:
            logging.debug("removing %s...", dirname)
            shutil.rmtree(dirname)

    def read_mast_output(self, mast_output, genes):
        """XML MAST output"""
        return mast_formats.from_xml_text(mast_output, genes)


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
    with tempfile.NamedTemporaryFile(mode='w+', prefix='memebg',
                                     delete=False) as outfile:
        filename = outfile.name
        #logging.info("make background file '%s'", filename)
        outfile.write("# %s order Markov background model\n" %
                      util.order2string(len(bgmodel) - 1))
        for order_row in bgmodel:
            for seq, frequency in order_row.items():
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


USER_TEST_FASTA_PATH = 'cmonkey/default_config/fasta_test.fa'

def is_meme_version_supported(version):
    if version is not None:
        major, minor, patch = map(int, version.split('.'))
        if major == 4:
            return True
        else:
            return False
    return True

def check_meme_version():
    logging.info('checking MEME...')
    try:
        test_fasta = resource_filename(Requirement.parse("cmonkey2"), USER_TEST_FASTA_PATH)
    except DistributionNotFound:
        test_fasta = USER_TEST_FASTA_PATH

    try:
        command = ['meme', '-nostatus', '-text', test_fasta]
        output = subprocess.check_output(command).decode('utf-8').split('\n')
        for line in output:
            if line.startswith('MEME version'):
                return line.split(' ')[2]
    except OSError:
        logging.error("MEME does not exist in your PATH, please either install or check your PATH variable")
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

def write_pssm(session, outfile, motif_info_id, evalue, num_sites):
    """writes a single PSSM to the given file"""
    outfile.write('\nMOTIF %d\n' % motif_info_id)
    outfile.write('BL   MOTIF %s width=0 seqs=0\n' % motif_info_id)

    pssm_rows = [(row.a, row.c, row.g, row.t)
        for row in session.query(cm2db.MotifPSSMRow).filter(cm2db.MotifPSSMRow.motif_info_id == motif_info_id)]

    outfile.write('letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e\n' % (len(pssm_rows), num_sites, evalue))
    for a, c, g, t in pssm_rows:
        outfile.write('%5.3f %5.3f %5.3f %5.3f\n' % (a, c, g, t))


def write_motifs2meme(session, filepath):
    """Write the motifs to a MEME file and returns True if successful
    Currently, this only works if there is global background data in the database
    """
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
    """a wrapper around the tomtom script"""
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
