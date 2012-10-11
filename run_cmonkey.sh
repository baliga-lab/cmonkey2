#!/bin/bash

PYTHONPATH=`pwd`/cmonkey:$PYTHONPATH python cmonkey/cmonkey.py $@

