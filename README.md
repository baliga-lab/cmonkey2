![cMonkey2 Logo](https://github.com/baliga-lab/cmonkey2/blob/master/graphics/cmonkey2_logo_80px.png "cMonkey2 Logo")

## cMonkey<sub>2</sub> - Python port of the [cMonkey biclustering algorithm](http://cmonkey.systemsbiology.net)

### Description

This is the Python implementation of the cMonkey algorithm based on the original R implementation by David J. Reiss, Institute for Systems Biology.

### Documentation

A complete set of documentation for installation and running of cMonkey is on the [wiki](https://github.com/baliga-lab/cmonkey2/wiki). There are also [developer](https://groups.google.com/d/forum/cmonkey-dev) and [user](https://groups.google.com/d/forum/cmonkey-users) discussion groups. 

### Contact

Please report all bugs or other issues using the [issue tracker](https://github.com/baliga-lab/cmonkey2/issues). Please direct any and all questions to either the [developer](https://groups.google.com/d/forum/cmonkey-dev) or [user](https://groups.google.com/d/forum/cmonkey-users) discussion groups. 

### Installation

The recommended way is to install cmonkey2 through pip

```
pip install cmonkey2
```

This will install the tools cmonkey2 and cm2view into your python environment. Please note that
you will have to install MEME manually from http://meme-suite.org/

### Running cmonkey2

The simplest way to run the tool (if all data available in RSAT and STRING):

```
$ cmonkey2 --organism <organism-code> <tab separated file of gene expressions>
```

To display available options:
```
bin/cmonkey2.sh --help
```

To run the example organism:
```
bin/cmonkey2.sh --organism hal --rsat_base_url http://networks.systemsbiology.net/rsat example_data/hal/halo_ratios5.tsv
```

### Using directly from the source repository

Below are the instructions to use cmonkey2 directly in the source repository

### System requirements

cMonkey<sub>2</sub> has been tested and runs on all tested recent versions of Linux (including debian-based [Ubuntu, Mint, Debian] and RPM-based [CentOS, Fedora]) and recent versions of Mac OS X. Additional dependencies include:

* Developed and tested with Python 2.7.x and Python 3.x
* scipy >= 0.9.0
* numpy >= 1.6.0
* biopython >= 1.63
* BeautifulSoup >= 4
* R >= 2.14.1
* rpy2 >= 2.2.1
* MEME 4.3.0 or >= 4.8.1
* csh (for running MEME)
* pandas
* sqlalchemy and sqlalchemy-utils
* svgwrite

for the human setup, Weeder 1.4.2 is needed

for running the unit tests (optional):

* python-xmlrunner

for running the interactive monitoring and visualization web application (optional):

* CherryPy 3
* Jinja2
* python-routes


### Running the Unit Tests

    bin/run_tests.sh

### Running cmonkey2

In general, you should be able to run cmonkey2 on microbial gene
expression ratios with

    bin/cmonkey2.sh --organism <organism-code> <tab separated file of gene expressions>

The file can be either in your file system or a web URL.

After the program was started, a log file will be written in cmonkey.log. You
can see all available options with

    bin/cmonkey2.sh --help

### Test Run with Halobacterium Salinarum

There is a startup script for cMonkey to run the current integrated
system

    bin/cmonkey2.sh --organism hal example_data/hal/halo_ratios5.tsv

### Start the python based monitoring application

    bin/cm2view.sh [--out [output directory]]

### Another way is to run Halobacterium is specify the RSAT database

    bin/cmonkey2.sh --organism hal --rsat_organism Halobacterium_NRC_1_uid57769 --rsat_base_url http://pedagogix-tagc.univ-mrs.fr/rsat --rsat_features gene --nooperons --use_BSCM example_data/hal/halo_ratios5.tsv


### Running cMonkey on Human

To run cMonkey on human data, run the following code with your own `<ratios.tsv>` file

    bin/cmonkey2.sh --organism hsa --string <stringFile> --rsat_organism Homo_sapiens_GRCh37 --rsat_URL http://rsat.sb-roscoff.fr/ --rsat_features protein_coding --nooperons <ratios.tsv>

#### More details for running cMonkey on human data

Running cMonkey on Human data is somewhat difficult because neither the string database nor the RSAT database has human data cleanly entered.  Here are the steps for a sucessful python cMonkey run on human

1.  Make a gene interaction file.  The example data file mentioned above was generated from Biogrid around 10/6/14.
2.  Find an RSAT mirror that has .raw chromose files and feature files.  In the above example, we use Homo\_sapiens\_ensembl\_74\_GRCh37 from the main RSAT database.  To annotate these we use 'protein\_coding.tab' and 'protein\_coding\_names.tab'.  In principal, other annotation files such as 'processed\_transcript' would work just as well.
3.  Adjust the upstream region searched, and perhaps modify the code to search for know TF and miRNA motifs rather than de-novo motifs.  NOTE: Modiyfing the motif search step is non-trivial.



### Package maintainers

#### General

The distribution is built using setuptools and wheel format

  - setup.py contains all information needed to build the distribution
    increase the version number before making a distribution
  - record user-relevant changes in CHANGELOG.rst

#### Build distribution

python3 setup.py sdist bdist_wheel

#### Uploading to PyPI

twine upload -r pypi dist/cmonkey2-<version>*
