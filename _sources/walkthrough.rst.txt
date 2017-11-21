Walkthrough of a Complete cMonkey\ :sub:`2` Run
===============================================

Setting up for a cMonkey\ :sub:`2` run
--------------------------------------

  * Input required The only required input are

    1. The gene expression matrix file (tab-separated values file)
    2. The KEGG organism code (e.g., eco for E. coli; hal for H. salinarum)

  * Input automatically downloaded Given the KEGG organism code, cMonkey\ :sub:`2` downloads the following data from various online databases and caches them locally:

    1. Genome sequences and annotations (currently from RSAT)
    2. Operon membership predictions (for prokaryotes)
    3. STRING functional gene associations

  * Input optional Custom files may be use as additional input (or to replace downloaded files):

    1. Network
    2. Genome
    3. Annotation

Starting the run
----------------

cMonkey\ :sub:`2` is typically started by running the front-end cmonkey2 with additional parameters.

Example:

.. highlight:: none

::

   cmonkey2 --organism hal --rsat_base_url http://networks.systemsbiology.net/rsat ratios.tsv

The following phases will be executed:

  1. Initialization
  2. Scoring
  3. Post-processing

Initialization steps
--------------------

Before cMonkey\ :sub:`2` runs the scoring iterations, it ensures that the required data is in place via the following steps:

  1. normalize the ratio matrix
  2. retrieve organism information and data from RSAT
  3. build gene synonym lookup table
  4. retrieve gene locations and parse promoter sequences
  5. download and normalize operon network from Microbes Online (if necessary)
  6. download and normalize STRING network (if necessary)
  7. initialize the starting biclusters

What is happening at each iteration
-----------------------------------

These are the essential steps in a cMonkey\ :sub:`2` iteration:

  1. run row scoring functions
  2. run column scoring function
  3. combine the scores via the combiner
  4. use scores to update cluster membership

Each scoring function can be individually configured to follow a specific schedule and the weights which define how much its results contribute to the final combined score can be updated following that schedule as well.

It should be noted that motif scoring has two schedules, one for running the MEME suite pipeline (typically every 100 iterations) and another one to apply the results of the MEME suite pipeline step to the current clusters memberships.

Inside a scoring module
-----------------------

The building blocks of a cMonkey\ :sub:`2` run are scoring modules. A scoring module is a Python class that implements the methods

.. code-block:: python

   def compute(self, iteration_result, reference_matrix)
   def compute_force(self, iteration_result, reference_matrix)

The only difference between these do functions is that ``compute_force()`` always performs a computation, while ``compute()`` only runs in the iterations it is scheduled for.

The result of these calls is a DataMatrix object, where the columns represent the clusters, and the rows represent the genes or conditions, depending on whether the scoring function is a row or a column scoring function. Each cell in the matrix contains the corresponding score for a gene/condition in a specific cluster.

cMonkey\ :sub:`2` comes with a number of built-in scoring modules, which are part of the standard scoring setup. It is possible to create user-defined scoring modules; in this case it is recommended to inherit from the ScoringFunctionBase class, which provides a lot of useful functionality.

Monitoring progress
-------------------

Throughout the computation, cMonkey\ :sub:`2` writes out its progress into its result database. Users can directly query the data contained in the database using regular data extraction tools. cMonkey\ :sub:`2` also comes with a web browser-based monitoring application that can be used to view a graphical representation of the run.

Post-processing steps
---------------------

Results data format
-------------------

Please see details about the result database schema here

Viewing results
---------------

Please see details about how to view results here

Plugging in cMonkey Results for downstream analysis
---------------------------------------------------
