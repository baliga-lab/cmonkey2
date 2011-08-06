"""meme.py - cMonkey meme suite integratino

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import subprocess


def dust(fasta_file_path):
    """runs the dust command on the specified FASTA file and
    returns a list of sequences"""
    output = subprocess.check_output(['dust', fasta_file_path])
    return output
