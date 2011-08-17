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


class MemeSuite:
    """Regard the meme suite as a unit of tools. This helps
    us capturing things like versions, global settings and data
    passing"""
    def __init__(self, max_width=24, use_revcomp=True):
        """Create MemeSuite instance"""
        self.__max_width = max_width
        self.__use_revcomp = use_revcomp

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
                logging.info("DUST input written to: %s", dust_input.name)
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
        return process_with_dust(seqs_for_dust)

    def run_meme(self, input_seqs, all_seqs):
        """Runs the meme tool. input_seqs is a dictionary of
        (feature_id : (location, sequence)) that are to be provided as meme
        input, all_seqs is a dictionary that provides all sequences used
        in the cMonkey run, which will be used to compute background
        distribution"""
        def add_if_unique(meme_input_seqs, seq):
            """add the sequence to the list only if it does not exist"""
            if seq not in meme_input_seqs:
                meme_input_seqs.append(seq)

        def make_seqs(seqs):
            """prepare the input sequences for feeding into meme.
            This means only taking the unique sequences and rever"""
            meme_input_seqs = []
            for locseq in seqs.values():
                seq = locseq[1]
                add_if_unique(meme_input_seqs, seq)
                if self.__use_revcomp:
                    add_if_unique(meme_input_seqs, st.revcomp(seq))
            return meme_input_seqs

        def background_seqs():
            """return all sequences to be used for background calculation"""
            return {feature_id: all_seqs[feature_id]
                    for feature_id in all_seqs if feature_id not in input_seqs}

        def make_background_file():
            """create a meme background file and returns its name"""
            bgseqs = background_seqs()
            bgfilename = None
            bgmodel = st.markov_background(make_seqs(bgseqs), 3)
            with tempfile.NamedTemporaryFile(prefix='memebg',
                                             delete=False) as outfile:
                bgfilename = outfile.name
                outfile.write("# %s order Markov background model\n" %
                              util.order2string(len(bgmodel) - 1))
                for order_row in bgmodel:
                    for seq, frequency in order_row.items():
                        outfile.write('%s %10s\n' %
                                      (seq, str(round(frequency, 8))))
            return bgfilename

        logging.info("run_meme() - # seqs = %d", len(input_seqs))
        bgfilename = make_background_file()
        logging.info("created background file in %s", bgfilename)

    def dust(self, fasta_file_path):
        """runs the dust command on the specified FASTA file and
        returns a list of sequences"""
        output = subprocess.check_output(['dust', fasta_file_path])
        return output

    def meme(self, input_file_path, background_file_path, psp_file_path,
             num_motifs):
        """runs the meme command on the specified input file, background file
        and positional priors file"""
        command = ("meme %s -bfile %s -psp %s -time 600 -dna -revcomp " +
                   "-maxsize 9999999 -nmotifs %d -evt 1e9 -minw 6 -maxw %d " +
                   "-mod zoops -nostatus -text") % (input_file_path,
                                                    background_file_path,
                                                    psp_file_path,
                                                    num_motifs,
                                                    self.__max_width)

    def mast(self, meme_outfile_path, database_file_path,
             background_file_path):
        """runs the mast command"""
        command = ("mast %s -d %s -bfile %s -nostatus -stdout -text -brief " +
                   "-ev 99999 -mev 99999 -mt 0.99 -seqp " +
                   "-remcorr") % (meme_outfile_path, database_file_path,
                                  background_file_path)
