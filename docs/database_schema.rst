Result Database Schema
======================

Introduction
------------

cmonkey-python writes the results of its computation to an SQLite database. This choice was made, because SQLite is a free, open source and portable data store which is available on many systems and has programming interfaces to a large number of programming languages. Another important aspect is that the entire database is stored in a single file, which can be easily copied, archived and analyzed. In this section the database structure and its function is explained in further detail.

Tables
------

**Note 1:** *The tables ending in _stats are only used in the cluster_viewer application and are subject to change.*

**Note 2:** *SQLite is different from other RDBMS in that each table has an implicit column rowid that acts like an auto incremented integer valued primary key. It is normally not shown in the frontend, but we will add it here for clarity*

run_infos
~~~~~~~~~

.. code-block:: sql

    rowid          int
    start_time     timestamp
    finish_time    timestamp
    num_iterations int
    last_iteration int
    organism       text
    species        text
    ncbi_code      int
    num_rows       int
    num_columns    int
    num_clusters   int
    git_sha        text

This table represents the current information about a cmonkey run and only stores a single entry that is continuously updated until a run is finished.

row_names, column_names
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: sql

    rowid     int
    order_num int
    name      text

These two tables are structurally identical. They reflect the structure of the input gene expression matrix, to preserve the order of the rows and columns, their order is stored as well.

row_members, column_members
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: sql

    rowid     int
    iteration int
    cluster   int
    order_num int

These tables contain the row and column members for each iteration and cluster. The element order_num references an order_num in its respective row_names/column_names table.

cluster_stats
~~~~~~~~~~~~~

.. code-block:: sql

    rowid     int
    iteration int
    cluster   int
    num_rows  int
    num_cols  int
    residual  decimal

Stores the residual values, number of rows and columns for each iteration and cluster.

motif_infos
~~~~~~~~~~~

.. code-block:: sql

    rowid     int
    iteration int
    cluster   int
    seqtype   text
    motif_num int
    evalue    decimal

Basic information about a motif that cmonkey thinks is associated with a specific cluster.

meme_motif_sites
~~~~~~~~~~~~~~~~

.. code-block:: sql

    rowid         int
    motif_info_id int  /* references motif_infos.rowid */
    seq_name      int
    reverse       boolean
    start         int
    pvalue        decimal
    flank_left    text
    seq           text
    flank_right   text

Detailed positional MEME information for a motif.

motif_annotations
~~~~~~~~~~~~~~~~~

.. code-block:: sql

    rowid         int
    motif_info_id int  /* references motif_infos.rowid */
    iteration     int
    gene_num      int
    position      int
    reverse       boolean
    pvalue        decimal

Positional information for a motif that was obtained from MAST.

motif_pssm_rows
~~~~~~~~~~~~~~~

.. code-block:: sql

    rowid         int
    motif_info_id int  /* references motif_infos.rowid */
    iteration     int
    row           int
    a             decimal
    c             decimal
    g             decimal
    t             decimal

Rows of the PSSM for a motif.

global_background
~~~~~~~~~~~~~~~~~

.. code-block:: sql

    rowid       int
    subsequence text
    pvalue      decimal

If the run uses a global background file, this table stores the entries that were generated.

statstypes
~~~~~~~~~~

.. code-block:: sql

    rowid      int
    category   text
    name       text

iteration_stats
~~~~~~~~~~~~~~~

.. code-block:: sql

    rowid      int
    statstype  int
    iteration  int
    score      decimal
