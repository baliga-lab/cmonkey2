cMonkey Principles
==================

In order to use cMonkey for a computation, a user needs to setup a **cMonkey run**. A run is a set of specific input data and configuration parameters and a set of scoring algorithms that are activated at certain iterations and combined using user-specified weights.

cMonkey Python supports the user by providing the CMonkeyRun class which contains a framework for executing a run. The user can override certain configuration parameters or its methods to customize the run.

Core algorithm
--------------

At the heart of the cMonkey algorithm is a scheduler that executes specific scoring modules at user-defined iterations and with a user-defined scoring weight. There are two distinct classes of scoring modules:

  * Row scoring modules are based on the input matrix's rows (or genes). The user may define an arbitrary number of row scoring modules in a run.
  * Column scoring modules are based on the input matrix's columns (or conditions). Currently, only a single column scoring function per run is supported

Anatomy of a scoring module
---------------------------

A cMonkey scoring module is a class that implements the method

.. code-block:: python

   def compute(self, iteration_result, reference_matrix=None)

which returns a DataMatrix object with

  * ``num_genes`` rows and ``num_clusters`` columns for a row scoring module where each cell contains the score for each gene-cluster combination
  * ``num_conditions`` rows and ``num_clusters`` columns for a column scoring module where each cell contains the score for each condition-cluster combination

System defined scoring modules
------------------------------

Row scoring modules
~~~~~~~~~~~~~~~~~~~

  * ``microarray.RowScoringFunction``: This is the standard row scoring function
  * ``network.ScoringFunction``: Uses a list of weighted networks to calculate scores
  * ``motif.MemeScoringFunction``: Calculates scores based on an external tool pipeline of meme and mast
  * ``motif.WeederScoringFunction``: Calculates scores based on an external tool pipeline of weeder and mast
  * ``scoring.ScoringFunctionCombiner``: A function module that represents an aggregate of one or more functions that are weighted and combined.

It is recommended to extend the base class ``scoring.ScoringFunctionBase``, to automatically inherit access to some useful and commonly used functionality.

Column scoring modules
~~~~~~~~~~~~~~~~~~~~~~

  * ``scoring.ColumnScoringFunction``: This is the default condition scoring function

How to implement your own scoring module
----------------------------------------

The simplest way to implement your own scoring module is to derive from ``scoring.ScoringFunctionBase``. This base class provides all necessary infrastructure for scheduling and basic configuration.
