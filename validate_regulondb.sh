#!/bin/bash

PYTHONPATH=`pwd`/cmonkey:$PYTHONPATH python cmonkey/tomtom_verify.py $@
tomtom -verbosity 2 -q-thresh 0.5 -dist ed -min-overlap 4 -text -query-pseudo 0 -target-pseudo 0 -target regulondb_pssms.txt -query meme_pssms.txt > tomtom_results.txt


