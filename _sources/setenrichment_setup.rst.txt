Set Enrichment Setup
====================

How to setup set enrichment
---------------------------

Set Enrichment is an optional scoring mechanism where the user can provide input files containing genes that are grouped in sets.

Input file Formats
~~~~~~~~~~~~~~~~~~

Supported set file formats are either JSON or CSV.

The JSON format is as follows:

.. highlight:: none

::

   {
       "set1": ["<gene1>", "<gene2>", ...],
       "set2": ["<gene1>", "<gene2>", ...],
       ...
   }

The CSV format is specified like this:

.. highlight:: none

::

   set1,<gene1>;<gene2>;...
   set2,<gene1>;<gene2>;...

Configuration
~~~~~~~~~~~~~
Users who want to use the set enrichment function need to specify it in a custom pipeline file

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
                   { "id": "SetEnrichment",
                       "function": { "module": "cmonkey.set_enrichment", "class": "ScoringFunction" }
                   },
                   ...
               ]
           }
       },
       ...
   }

In addition, the user needs to provide the Enrichment settings in an .ini file:

.. highlight:: none

::

   [SetEnrichment]
   schedule = 1,7
   scaling_rvec=seq(1e-5, 0.5, length=num_iterations*3/4)
   set_types = settype1

   [SetEnrichment-settype1]
   set_file = example_data/eco/eco-merged-set.json
   weight = 1.0

As shown, the ``SetEnrichment`` section specifies the schedule and the scaling within the row scoring functions as well as a single set type. The user can specify any number of set types, by providing a list of names in the set_types setting, separated by comma. Set type specific settings are then made in the section ``"SetEnrichment-"``, currently these are the set file and the weight that indicates the percentage that scoring on a set type contributes to set enrichment scoring.
