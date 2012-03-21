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

