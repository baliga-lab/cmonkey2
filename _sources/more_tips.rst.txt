More Tips for Running, Debugging, Developing
============================================

Interactive exploration of results
----------------------------------

Via cm2view
~~~~~~~~~~~

Run statistics, bicluster statistics, views of individual biclusters, and a network visualization are available using the interactive cMonkey2 viewer.

Via IPython notebook
~~~~~~~~~~~~~~~~~~~~

The IPython interactive Python console enables interactive exploration of cMonkey2 results via an interactive notebook-style interface. This is further facilitated by a cmonkeyobj python module (still in development).

For an example of what is possible, check out this example IPython notebook.

Note that cMonkey2 may also be run in an IPython notebook. While this can take a long time, it will also (eventually) enable cMonkey2 input and output to be captured and recalled. We have a preliminary example of this capability here.

Debugging
---------

Interactive Debugging
~~~~~~~~~~~~~~~~~~~~~

The IPython interactive Python console also provides a variety features that do not exist in the default Python console. It can be useful for interactive work, debugging, testing, or developing with cmonkey-python. It provides an interactive environment that is similar to the R console for cMonkeyR.

If either python or ipython are called with the -i switch, cmonkey will remain in the interpreter after running, which can facilitate debugging. This can be combined with cmonkey's ``--interactive`` switch, e.g.

.. highlight:: none

::

   ipython -i bin/cmonkey2 -- --organism hal example_data/hal/halo_ratios5.tsv --interactive.

This will initialize the run, but stop before scoring, making the object cmonkey_run available for interactive access.

Data exchange
-------------

  * A utility for importing JSON output into R: http://cran.r-project.org/web/packages/RJSONIO/index.html
