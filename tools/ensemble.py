#!/usr/bin/env python
"""ensemble.py
This is a tool to prepare a directory for ensemble runs. Given an organism code
and the ratios, matrix, split up the matrix into sub matrices and create
SGE qsub shell scripts
"""
import argparse
import os
import datamatrix as dm

DESCRIPTION = """ensemble.py - prepare cluster runs"""

QSUB_TEMPLATE = """#!/bin/bash

export LD_LIBRARY_PATH=/tools/lib:/tools/R-3.0.3/lib64/R/lib
export PATH=/tools/bin:${PATH}

#$ -S /bin/bash
#$ -m be
#$ -q baliga
#$ -P Baliga
#$ -M wwu@systemsbiology.org
#$ -cwd
#$ -pe serial 8
#$ -l mem_free=32G

python cmonkey.py --organism %s --ratios %s --out %s"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True)
    parser.add_argument('--ratios', required=True)
    parser.add_argument('--targetdir', required=True)
    parser.add_argument('--numfiles', default=4)
    parser.add_argument('--numcols', default=8)

    args = parser.parse_args()
    dm.prepare_ensemble_matrix(args.ratios, args.targetdir, args.numfiles, args.numcols)
    for i in range(1, args.numfiles + 1):
        with open(os.path.join(args.targetdir, "%s-%03d.sh" % (args.organism, i)), 'w') as outfile:
            outfile.write(QSUB_TEMPLATE % (args.organism,
                                           os.path.join(args.targetdir, "ratios-%03d.tsv.gz" % i),
                                           "%s-out-%03d" % (args.organism, i)))
