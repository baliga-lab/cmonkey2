# Inferelator Integration package

This directory contains tools to integrate with the Inferelator version at https://github.com/dreiss-isb/cMonkeyNwInf

The script to use is run_inf.R. It can be used in 3 ways:


  1. provide a list of transcription factors and a cmonkey-python
     output directory.
  Example:
```
        ./run_inf.R --tfsfile <tfsfile> --resultdir <resultdir> --outfile <outfile>
```

  2. provide a list of transcription factors, a gzip compressed ratios matrix and
     a cluster stack in JSON format
  Example:
```
./run_inf.R --tfsfile <tfsfile> --json <jsonfile> --ratios <ratios> --outfile <outfile>
```

The JSON format for cluster stacks looks like this:

```
[
  {
    "nrows": <# of row members>,
    "ncols": <# of column members>,
    "rows": [ <list of row members> ],
    "cols": [ <list of column members> ],
    "k": <cluster number>,
    "resid": <cluster residual>
  },
  ...
]    
```

  3. provide an R environment as an RData file that contains the
     following parameters:
     * ratios: the gene expression matrix
     * predictors: the list of transcription factors
     * e: an environment with the following data
     * clusterStack: cMonkey standard cluster stack format
```
Example: ./run_inf.R --rdata <rdatafile> --outfile <outfile>
```

