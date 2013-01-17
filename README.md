## cMonkey Python - Python port of the cMonkey biclustering algorithm

### Description

This is the Python implementation of the cMonkey algorithm based on the original R implementation by David Reiss, Institute for Systems Biology

### Status

Testing

### System requirements

* Developed and tested with Python 2.7.2
* scipy >= 0.9.0
* numpy >= 1.6.0
* MySQLdb >= 1.2.3
* BeautifulSoup >= 3.2.0
* R >= 2.14.1
* rpy2 >= 2.2.1
* MEME 4.3.0
* csh (for running MEME)
for the human setup, Weeder 1.4.2 is needed

### Running the Unit Tests

    ./run_tests.sh


### Running cMonkey

In general, you should be able to run cmonkey-python on microbial gene
expressions with

    ./run_cmonkey.sh --organism <organism-code> --ratios <tab separated file of gene expressions>

The file can be either in your file system or a web URL.

### Test Run with Halobacterium Salinarum

There is a startup script for cMonkey to run the current integrated
system

    ./run_cmonkey.sh --organism hal --ratios halo_ratios5.tsv

