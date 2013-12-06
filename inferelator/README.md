# Inferelator Integration package

This directory contains tools to integrate with the Inferelator version at https://github.com/dreiss-isb/cMonkeyNwInf

The script to use is run_inf.R. Please look at the header comments in this file for usage
details.

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