Input File Formats
==================

cMonkey uses a number of input files to support its built in scoring functions.

**Note:** cMonkey-Python uses a cache directory for data files that are downloaded from the web. You can use this mechanism for specifying your own data, as long as it is in the same format as specified below.

The path is by default ``cache``, you can also specify a different location by providing the ``--cachedir`` option on the command line

Gene expressions (mandatory)
----------------------------

The most basic and mandatory input is gene expression data. This should be specified as a tab-delimited text file, with the following specifications for a matrix with n genes and m conditions:

  * first row: 1 dummy label + m condition titles
  * n rows where each row starts with a gene name and m expression ratio values

The expression file can either be uncompressed or in gzip format (in this case, it should have a .gz suffix).

**Example:**

.. highlight:: none

::

  X<TAB>Cond 1<TAB>Cond 2
  Gene1<TAB>1.213<TAB>-1.412
  ...

RSAT information
----------------

The RSAT database is central to cMonkey's automatic retrieval of organism information. It is used for several purposes:

  * organism classification (eukaryote/prokaryote, NCBI taxonomy mapping)
  * gene synonyms
  * genomic information

RSAT provides organism data on a number of sites, where there is usually a web directory under a URL that looks like ``<BASE_URL>/data/genomes/``

This directory contains many organism directories. The user can specify this RSAT organism explicitly via the ``--rsat_organism`` parameter on the command line, otherwise, it will attempt to derive the name through the mandatory KEGG organism code. In addition the base url of the RSAT site to be searched can be defined with the ``--rsat_base_url`` parameter.

cMonkey\ :sub:`2` expects the following files in an organism directory:

RSAT organism file ("organism.tab")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The contents of this file is used to look for the occurrence of the word "Eukaryota" and derive the NCBI taxonomy id. If the word "Eukaryota" is found, the organism will be marked as a eukaryote, otherwise, it will be assumed as prokaryote.

**Example:**

.. highlight:: none

::

   64091<TAB>Archaea; Euryarchaeota; Halobacteria; Halobacteriales; Halobacteriaceae; Halobacterium

RSAT feature name files ("feature_names.tab")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These files are used by cMonkey to create the synonym table for alternative gene names. Each line has the format

.. highlight:: none

::

   <Accession ID><TAB><name><primary|alternate>

**Example:**

.. highlight:: none

::

  NP_045946.1<TAB>VNG7001<TAB>primary
  NP_045946.1<TAB>1446803<TAB>alternate
  NP_045946.1<TAB>10803548<TAB>alternate
  NP_045946.1<TAB>NP_045946.1<TAB>alternate
  ...

RSAT feature files ("feature.tab")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each line in a feature file contains a gene location for an organism. The important information is are the id, name, contig, strand, start and end position. For each contig listed in this file, there must exist a corresponding contig file.

**Example:**

.. highlight:: none

::

   -- id<TAB>type<TAB>name<TAB>contig<TAB>start<TAB>end<TAB>strand<TAB>description<TAB>chrom_pos<TAB>organism<TAB>GeneID
   NP_045946.1<TAB>CDS<TAB>VNG7001<TAB>NC_001869.1<TAB>363<TAB>812<TAB>R<TAB>hypothetical protein<TAB>complement(363..812)<TAB>Halobacterium sp. NRC-1<TAB>1446803
   NP_045947.1<TAB>CDS<TAB>VNG7002<TAB>NC_001869.1<TAB>834<TAB>1172<TAB>R<TAB>hypothetical protein<TAB>complement(834..1172)<TAB>Halobacterium sp. NRC-1<TAB>1446804
   ...

RSAT contig files
~~~~~~~~~~~~~~~~~

These files contain the raw genomic sequence in lower case for a specific contig/chromosome that is referenced in the RSAT features file.

**Example (file name "NC002607.1.tab"):**

.. highlight:: none

::

   ttgacccactgaatcacgtctgaccgcgcgtacgcggtcacttgcggtgccgttttctttgttaccgacgaccgaccagcgacagccaccgcgcgctcactgccaccaaaagagtcatatcacagccgaccagtttctggaacgttcccgatactggaacggtcctaatgcagtatcccaccctccttccatcgacgccagtcgaatcacgccgccagccaccgtccgccagccggccagaataccgatgactcggcggtctcgtgtcggtgccggcctcgcagccattgtactggccctggccgcagtgtcggctgccgctcc

RSAT mockup directories
~~~~~~~~~~~~~~~~~~~~~~~

For a number of organisms, RSAT does not provide the necessary files for a cmonkey run. In this case, the user can provide their own local directory that mirrors the structure, so it contains the files described above and specify the directory using --rsat_dir <your_rsat_directory>

**Example:**

.. highlight:: none

::

   organism.tab
   feature_names.tab
   feature.tab
   <contig1>.tab
   <contig2>.tab
   ...

STRING protein-protein interactions
-----------------------------------

STRING is a database of known and predicted protein-protein interactions. cMonkey uses these interactions to improve its clustering results. In order to do so, it builds a network using the interactions for the genes for the current organism. This network can then be used in the network scoring component of cMonkey.

STRING is an enormous database and so it is a good idea to prepare the input to provide only the necessary data for cMonkey's network scoring algorithm, namely the names of the genes and their score.

STRING files are tab-delimited files containing entries of the form

.. highlight:: none

::

   Gene1<TAB>Gene2<TAB>Normalized Score

cmonkey-python provides the utility ``extract_string_links.sh`` to write the interactions for an organism, given the KEGG code and a database file (either gzip'ed or plain).

Protein-protein interaction network files have a very simple structure, each line in the file defines a a network edge with a weight, and each file is compressed with gzip, so the naming scheme is <NCBI code>.gz.

**Example (file name "64091.gz"):**

.. highlight:: none

::

  VNG0001H<TAB>VNG0002G<TAB>865
  VNG0001H<TAB>VNG0003C<TAB>802
  VNG0001H<TAB>VNG0005H<TAB>561
  ...

Microbes Online operon files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These files are typically provided by Microbes Online and define operon relationships of genes. These are only used for prokaryotic organisms. The name scheme that cmonkey uses for such files is gnc<NCBI code>.named. See example for the format of the lines in a file.

**Example (file name "gnc64091.named"):**

.. highlight:: none

::

   Gene1<TAB>Gene2<TAB>SysName1<TAB>SysName2<TAB>Name1<TAB>Name2<TAB>bOp<TAB>pOp<TAB>Sep<TAB>MOGScore<TAB>GOScore<TAB>COGSim<TAB>ExprSim
   68130<TAB>68131<TAB>VNG0001H<TAB>VNG0002G<TAB>VNG0001H<TAB>yvrO<TAB>TRUE<TAB>0.950<TAB>-3.000<TAB>0.000<TAB>NA<TAB>N<TAB>0.791
   68131<TAB>68132<TAB>VNG0002G<TAB>VNG0003C<TAB>yvrO<TAB>VNG0003C<TAB>TRUE<TAB>0.920<TAB>30.000<TAB>0.600<TAB>1.000<TAB>Y<TAB>0.525
   ...

