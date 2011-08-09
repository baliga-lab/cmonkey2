"""meme.py - cMonkey meme suite integration
These are functions and classes to conveniently call the needed
commands of the MEME suite in order to find motifs

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import subprocess


def dust(fasta_file_path):
    """runs the dust command on the specified FASTA file and
    returns a list of sequences"""
    output = subprocess.check_output(['dust', fasta_file_path])
    return output


def meme(input_file_path, background_file_path, psp_file_path, num_motifs):
    """runs the meme command on the specified input file, background file
    and positional priors file"""
    command = ("meme %s -bfile %s -psp %s -time 600 -dna -revcomp " +
               "-maxsize 9999999 -nmotifs %d -evt 1e9 -minw 6 -maxw 24" +
               "-mod zoops -nostatus -text") % (input_file_path,
                                                background_file_path,
                                                psp_file_path,
                                                num_motifs)


def mast(meme_outfile_path, database_file_path, background_file_path):
    command = ("mast %s -d %s -bfile %s -nostatus -stdout -text -brief " +
               "-ev 99999 -mev 99999 -mt 0.99 -seqp " +
               "-remcorr") % (meme_outfile_path, database_file_path,
                              background_file_path)
