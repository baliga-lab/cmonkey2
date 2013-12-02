#!/bin/bash
RATIOS="tmp/ecoli.tsv"
INDIR="ensemble_in"
N=2
K=8

PYTHONPATH=cmonkey python -c "import datamatrix as dm; dm.prepare_ensemble_matrix(\"$RATIOS\", \"$INDIR\", $N, $K)"

for f in ensemble_in/*.tsv.gz ; do
  filename="${f##*/}"
  filename=${filename%.tsv.gz}
  seq="${filename##ratios}"
  outdir="out$seq"
  echo "./cmonkey.py --organism eco --ratios $f --out $outdir --config example_data/eco/ecoli.ini"
  ./cmonkey.py --organism eco --ratios $f --out $outdir --config example_data/eco/ecoli.ini
done

