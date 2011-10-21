"""weeder.py - cMonkey weeder interface module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import subprocess as sp
import logging


LAUNCHER = 'weederlauncher'

class Weeder:
    def __init__(self):
        pass

    def run(self, fasta_file):
        command = [LAUNCHER, fasta_file, 'HS3P', 'small', 'T50']
        with open('weeder.log', 'w') as logfile:
            logging.info("running weeder on '%s'", fasta_file)
            weederproc = sp.Popen(command, stdout=logfile, stderr=sp.STDOUT)
            retcode = weederproc.wait()
            logging.info("Weeder finished, return code: %d", retcode)
