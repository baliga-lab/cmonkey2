#!/usr/bin/env python
"""ensemble.py
This is a tool to prepare a directory for ensemble runs. Given an organism code
and the ratios, matrix, split up the matrix into sub matrices and create
SGE qsub shell scripts

Example:

PYTHONPATH=cmonkey tools/ensemble.py --organism eco --ratios ecoli.tsv.gz --targetdir eco-ens-20140602 --numfiles 20 --mincols 50 --num_cores 1
"""
import argparse
import os
import datamatrix as dm

DESCRIPTION = """ensemble.py - prepare cluster runs"""

# Templates for Bourne Shell
QSUB_TEMPLATE_HEADER = """#!/bin/bash

export LD_LIBRARY_PATH=/tools/lib:/tools/R-3.0.3/lib64/R/lib
export PATH=/tools/bin:${PATH}
export BATCHNUM=`printf "%03d" $SGE_TASK_ID`
"""

QSUB_TEMPLATE = """#$ -S /bin/bash
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -t 1-%d
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=32G

python cmonkey_ensemble.py --organism %s --ratios %s --out %s --num_cores %d --ensemble_run_id $SGE_TASK_ID --minimize_io"""


# Templates for csh

QSUB_TEMPLATE_HEADER_CSH = """#!/bin/csh

setenv LD_LIBRARY_PATH /tools/lib:/tools/R-3.0.3/lib64/R/lib
setenv PATH /tools/bin:${PATH}
setenv BATCHNUM `printf "%03d" $SGE_TASK_ID`
"""

QSUB_TEMPLATE_CSH = """#$ -S /bin/csh
#$ -m be
#$ -q baliga
#$ -P Bal_%s
#$ -t 1-%d
#$ -M %s@systemsbiology.org
#$ -cwd
#$ -pe serial %d
#$ -l mem_free=32G

python cmonkey_ensemble.py --organism %s --ratios %s --out %s --num_cores %d --ensemble_run_id $SGE_TASK_ID"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--organism', required=True)
    parser.add_argument('--ratios', required=True)
    parser.add_argument('--targetdir', required=True)
    parser.add_argument('--numfiles', type=int, default=4)
    parser.add_argument('--mincols', type=int, default=8)
    parser.add_argument('--num_cores', type=int, default=1)
    parser.add_argument('--user', default=None)
    parser.add_argument('--csh', action='store_true')
    args = parser.parse_args()

    if args.csh:
        header = QSUB_TEMPLATE_HEADER_CSH
        template = QSUB_TEMPLATE_CSH
    else:
        header = QSUB_TEMPLATE_HEADER
        template = QSUB_TEMPLATE

    dm.prepare_ensemble_matrix(args.ratios, args.targetdir, args.numfiles,
                               args.mincols)
    with open(os.path.join(args.targetdir, "%s.sh" % args.organism), 'w') as outfile:
        if args.user is not None:
            login = args.user
        else:
            login = os.getlogin()

        outfile.write(header)
        outfile.write(template % (login, args.numfiles, login,
                                  args.num_cores,
                                  args.organism,
                                  os.path.join(args.targetdir, "ratios-$BATCHNUM.tsv.gz"),
                                  "%s-out-$BATCHNUM" % (args.organism),
                                  args.num_cores))
