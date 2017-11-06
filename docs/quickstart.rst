Quick Start
===========

After :doc:`installing <install>` the cmonkey2 package and ensuring that `MEME <http://meme-suite.org/>`_
is installed and in the search path, download the
`example data <http://networks.systemsbiology.net/downloads/cmonkey2-examples/>`_. This
is a simple gene expression matrix in the format that cmonkey2 expects.

.. highlight:: none

::

   cmonkey2 --organism hal --rsat_base_url http://networks.systemsbiology.net/rsat halo_ratios5.tsv


If everything is installed correctly, cmonkey2 will now take some time until
the run finishes.

You can see the progress of your run by running

.. highlight:: none

::

   cm2view

in the same directory as you started the run. This will start the web application
"cmonkey2 viewer", which can be accessed at the default address http://localhost:8080.
The browser should display something similar to this:

.. image:: _static/cm2view-screenshot-01.png

After cmonkey2 has finished the run, the results will be available in the ``out`` directory:

  * ``ratios.tsv.gz`` - the normalized input matrix
  * ``cmonkey_run.db`` - the results of the run as a sqlite3 database file
