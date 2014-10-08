## cMonkey Python - Python port of the cMonkey biclustering algorithm

### Description

This is the Python implementation of the cMonkey algorithm based on the original R implementation by David J. Reiss, Institute for Systems Biology.

### Documentation

A complete set of documentation for installation and running of cMonkey is on the [wiki](https://github.com/weiju/cmonkey-python/wiki). There are also [developer](https://groups.google.com/d/forum/cmonkey-dev) and [user](https://groups.google.com/d/forum/cmonkey-users) discussion groups. 

### System requirements

* Developed and tested with Python 2.7.2
* scipy >= 0.9.0
* numpy >= 1.6.0
* biopython >= 1.63
* MySQLdb >= 1.2.3
* BeautifulSoup >= 3.2.0
* R >= 2.14.1
* rpy2 >= 2.2.1
* MEME 4.3.0 or >= 4.8.1
* csh (for running MEME)
* Biopython >= 1.6.4 (earlier may work)
for the human setup, Weeder 1.4.2 is needed

for running the monitoring application (optional):

* CherryPy 3
* Jinja2
* python-routes

### Running the Unit Tests

    ./run_tests.sh

### Running cMonkey

In general, you should be able to run cmonkey-python on microbial gene
expressions with

    ./cmonkey.py --organism <organism-code> --ratios <tab separated file of gene expressions>

The file can be either in your file system or a web URL.

After the program was started, a log file will be written in cmonkey.log. You
can see all available options with

    ./cmonkey.py --help

### Test Run with Halobacterium Salinarum

There is a startup script for cMonkey to run the current integrated
system

    ./cmonkey.py --organism hal --ratios example_data/hal/halo_ratios5.tsv

### Start the python based monitoring application

    python cmviewer/main.py

### Running cMonkey on Human

To run cMonkey on human data, run the following code with your own `<ratios.tsv>` file

    ./cmonkey.py --organism hsa --ratios <ratios.tsv> --string ./example_data/hsa/BIOGRID-ALL-3.2.117.string.tsv.gz --rsat_organism Homo_sapiens_ensembl_74_GRCh37 --rsat_URL http://rsat.bigre.ulb.ac.be/rsat/ --rsat_features protein_coding --nooperons

#### More details for running cMonkey on human data

Running cMonkey on Human data is somewhat difficult because neither the string database nor the RSAT database has human data cleanly entered.  Here are the steps for a sucessful python cMonkey run on human

1.  Make a gene interaction file.  The example data file mentioned above was generated from Biogrid around 10/6/14.
2.  Find an RSAT mirror that has .raw chromose files and feature files.  In the above example, we use Homo\_sapiens\_ensembl\_74\_GRCh37 from the main RSAT database.  To annotate these we use 'protein\_coding.tab' and 'protein\_coding\_names.tab'.  In principal, other annotation files such as 'processed\_transcript' would work just as well.
3.  Adjust the upstream region searched, and perhaps modify the code to search for know TF and miRNA motifs rather than de-novo motifs.  NOTE: Modiyfing the motif search step is non-trivial.
