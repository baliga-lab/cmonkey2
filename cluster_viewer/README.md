## cMonkey Cluster Viewer

This is a web application for viewing the results of a cMonkey run.
It takes the following pieces of data:

- a JSON output directory
- a synonyms file and format
- a ratios file

and produces an HTML view of the current results

### System requirements

This application was written in Scala and the Play framework 2.0
(http://www.playframework.org)

In order to run it, it is necessary to install Play framework 2.0, set
the execution path and in the cluster_viewer directory enter

play
run

The configuration requires the following information:

- output directory containing the JSON files
- gene expression ratios
- synonyms

The user can provide this information in the project.conf file.

### Basic instructions for viewing the halobacteria test sequence

These instructions worked on Ubuntu 11.10, after running './run_cmonkey.sh hal halo_ratios5.tsv string_links_64091.tab'

1. Download the Play 2.0 framework (http://www.playframework.org) and unzip the play-2.0 directory and put it somewhere (e.g. /nethome/sdanziger/Programs/play-2.0).

2. Change into the cluster_viewer directory (e.g. cd /nethome/sdanziger/Work/Phospho/pMonkeyPython/cmonkey-python/cluster_viewer).

3. Edit 'project.conf' and change 'cmonkey.out.directory' to a directory that contains the *.json files and 'ratios.tsv' (e.g. cmonkey.out.directory=/nethome/sdanziger/Work/Phospho/pMonkeyPython/cmonkey-python/out).

4. Run Play 2.0 (e.g. /nethome/sdanziger/Programs/play-2.0/play)

5. Type "run"

6. Point a web browser towards "http://localhost:9000/"

