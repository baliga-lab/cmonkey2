Visualize Results in cMonkey2 Viewer
====================================

Description
-----------

You can view your current running or completed cMonkey2 run using cm2view. This is a web application that reads the output directory and generates a graphical representation from the data currently presiding there. As described in the (quickstart)[quickstart], cmviewer may be initialized by running

.. highlight:: none

::

    cm2view --out <your result directory>

The monitoring application will be started and can be viewed by opening a web browser to the address http://localhost:8080.

Views
-----

The views provided are:

  * Run statistics: Mean scores, cluster membership distributions etc.
  * Network view
  * List of clusters with their residuals and motif e-values
  * Cluster statistics: Genes and conditions in the cluster as well as motifs and their locations

Moreover, individual clusters in the cluster list may be clicked to visualize their data (including gene expression, motifs, motif positions, etc.).

Gaggle microformats for integration with FireGoose/ChromeGoose are included in each page.
