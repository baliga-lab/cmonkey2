#!/bin/bash

echo "Running unit tests..."
PYTHONPATH=`pwd`/cmonkey python test/all_tests.py $@
