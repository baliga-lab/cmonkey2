#!/bin/bash

PYTHONPATH=`dirname $0`/..
PYTHONPATH=$PYTHONPATH python `dirname $0`/addnwportal.py $@
