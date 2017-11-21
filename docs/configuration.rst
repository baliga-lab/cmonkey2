Configuration and Run Options
=============================

cMonkey\ :sub:`2` retrieves its configuration parameters through 2 mechanisms:

  1. a hierarchy of configuration files
  2. command line options

When the application is started, it will load its default configuration file ``default.ini``. It is recommended to leave the settings in this file unmodified.

Default configuration settings may be overridden by the user with a custom (set of) configuration file(s) and providing the ``--config <config file 1> <config file 2>...`` command-line option. Only the settings that should be changed need to be provided. Note that it is possible to specify multiple configuration files. In this case, the settings contained in the separate configuration files will be merged.

Some parameters can also be provided on the command line. In general, the settings in the user provided configuration file override the default parameters and command line parameters override configuration file parameters.

Command line options
--------------------

A current list of available command line options and their meaning in ``cmonkey.py`` can be retrieved by entering

.. highlight:: none

::

   cmonkey2 --help

Configuration (.ini) files
--------------------------

``cmonkey.py`` uses ``.ini`` files as a generic configuration mechanism for global and scoring specific settings. Global settings are defined in the ``[General]``, ``[Membership]`` and ``[Scoring]`` sections while each scoring function can implement its own set of specific configuration settings. The specific sections for each scoring functions are typically named after the name of the scoring functions.

Useful configuration constructs
-------------------------------

Scheduling
~~~~~~~~~~

Scheduling is a central part of cMonkey\ :sub:`2`. It defines in which iterations certain parts of the scoring algorithm will be active. For example, MEME is not run in every iteration, because its invocation is very costly. Instead, it is usually defined to run every 100 iterations or so. In order to define MEME to run every 100 iterations starting from iteration 300, we would define it like this:

.. highlight:: none

::

   [MEME]
   schedule=300:100

Or the user might only want to call MEME in iteration 42:

.. highlight:: none

::

   [MEME]
   schedule=42

We can define even more flexible schedules by combining a set of sub-schedules within a schedule by separating the sub-schedules with the '``:``' character.

.. highlight:: none

::

   [MEME]
   schedule=15:42:300,100

The above schedule says "run MEME in iteration 15, 42 and starting from 300, every 100 iterations". These principles can be applied to any setting that defines a schedule, e.g. network scoring, row scoring etc.

Scaling
~~~~~~~

The user can specify the strength of the contribution of each scoring module's scores to in each iteration. This can be done through the scaling setting. Currently two ways are allowed:

  1. ``scaling_const`` specifies a constant weight, given as a number in each iteration over the entire time of the run
  2. ``scaling_rvec`` allows the user to specify an R expression that can use the num_iterations variable to generate a vector of weights for each iteration

Scoring pipeline configuration (.json) files
--------------------------------------------

By default, cmonkey uses the default_pipeline.json file to setup a standard scoring pipeline. Users can choose to override the scoring functions used (e.g. for using user-defined scoring functions) by providing their own configuration files.

The general structure of a pipeline configuration file is as documented below:

.. highlight:: none

::

   {
     "row-scoring": {
       "id": "combiner",
       "function": { "module": "cmonkey.scoring", "class": "ScoringFunctionCombiner" },
       "args": {
         "functions": [
           { "id": "Rows",
             "function": { "module": "cmonkey.microarray", "class": "RowScoringFunction" }
           },
           ...
         ]
       }
      },
      "column-scoring": { "id": "Columns",
                          "function": { "module": "cmonkey.scoring",
                                        "class": "ColumnScoringFunction"} }
    }

The current version of the software is tested using a combiner function for row scoring and and a single scoring function for column scoring and this is the recommended setup. It is left to the users to pick the set of scoring functions most suitable to the problem and to provide this in their own .json files. With the --pipeline switch, the the path to the user-defined configuration file can be specified on the command line.

Tricky topics
-------------

The cluster is scored very good, but there are no motifs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cmonkey2 will not run MEME if the number of unique sequences for a cluster is less than a certain threshold value. This value is by default configured as

.. highlight:: none

::

   [Membership]
   ...
   min_cluster_rows_allowed = 3

If this hides good motifs that otherwise should be there, this number can be set to a lower value to ensure MEME is run on these cluster sequences.
