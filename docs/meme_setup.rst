MEME Setup
==========

MEME Suite notes
----------------

The tested and supported MEME versions are 4.3.0, 4.8.x, 4.9.0 and 4.10.x. cmonkey-python can automatically detect the installed version, provided it is in the executable search path. It is recommended to apply all available patches for the used MEME suite to ensure maximum stability. Please see the installation instructions for MEME for further details.

Installation script
-------------------

There is a small shell script that will download, compile and install MEME 4.10.0_2 (the latest as of this writing) in a local ``./progs`` directory, and make symlinks to the relevant tools in the current directory. This will enable local installation without requiring root. To run:

.. highlight:: none

::

   cmonkey2/tools/install_meme.sh

MEME 4.3.0
----------

To install this MEME version, the following procedure was used:

.. highlight:: none

::

  wget http://ebi.edu.au/ftp/software/MEME/4.3.0/meme_4.3.0.tar.gz
  tar -xvzf meme_4.3.0.tar.gz; cd meme_4.3.0
  ./configure --prefix=<your prefix> --enable-dependency-tracking --enable-opt --disable-shared --disable-fast-install --enable-serial --enable-build-libxml2 --enable-build-libxslt --enable-static --with-gnu-ld
  make
  make test
  sudo make install

MEME 4.8.x and later
--------------------

Starting from MEME versions 4.8.x, output data formats and switches are slightly different and the settings for building and installation are simpler:

.. highlight:: none

::

  tar -xvzf meme_<version>.tar.gz; cd meme_<version>
  ./configure --prefix=<your prefix>
  make
  make test
  sudo make install

Place ``<your prefix>/bin`` in your search path so cmonkey can find the MEME tools. You can of course, set the prefix to be anywhere else, as long as you make sure that the MEME tools can be found in the search path later.
