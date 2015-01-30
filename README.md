![cMonkey2 Logo](https://github.com/baliga-lab/cmonkey2/blob/master/graphics/cmonkey2_logo_80px.png "cMonkey2 Logo")

## cMonkey<sub>2</sub> - Python port of the [cMonkey biclustering algorithm](http://cmonkey.systemsbiology.net)

### Description

This is the Python implementation of the cMonkey algorithm based on the original R implementation by David J. Reiss, Institute for Systems Biology.

### Documentation

A complete set of documentation for installation and running of cMonkey is on the [wiki](https://github.com/baliga-lab/cmonkey2/wiki). There are also [developer](https://groups.google.com/d/forum/cmonkey-dev) and [user](https://groups.google.com/d/forum/cmonkey-users) discussion groups. 

### System requirements

cMonkey<sub>2</sub> has been tested and runs on all tested recent versions of Linux (including debian-based [Ubuntu, Mint, Debian] and RPM-based [CentOS, Fedora]) and recent versions of Mac OS X. Additional dependencies include:

* Developed and tested with Python 2.7.x
* scipy >= 0.9.0
* numpy >= 1.6.0
* biopython >= 1.63
* MySQLdb >= 1.2.3
* BeautifulSoup >= 3.2.0
* R >= 2.14.1
* rpy2 >= 2.2.1
* MEME 4.3.0 or >= 4.8.1
* csh (for running MEME)
for the human setup, Weeder 1.4.2 is needed

for running the unit tests (optional):

* python-xmlrunner 

for running the interactive monitoring and visualization web application (optional):

* CherryPy 3
* Jinja2
* python-routes

### Running the Unit Tests

    ./run_tests.sh

### Running cmonkey2

In general, you should be able to run cmonkey2 on microbial gene
expression ratios with

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
