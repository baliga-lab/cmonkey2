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


class MemeSuite:
    """Regard the meme suite as a unit of tools. This helps
    us capturing things like versions, global settings and data
    passing"""
    def __init__(self, max_width=24):
        """Create MemeSuite instance"""
        self.__max_width = max_width

    def remove_low_complexity(self, seqs):
        """send sequences through dust filter, send only those
        to dust that are larger than max_width"""
        def process_with_dust(seqs):
            """data conversion from and to dust tool"""
            dust_tmp_file = None
            with tempfile.NamedTemporaryFile(delete=False) as dust_input:
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
