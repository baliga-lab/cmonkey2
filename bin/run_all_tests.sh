#!/bin/bash

echo "Running unit tests..."
PYTHONPATH=`pwd` test/all_tests.py $@
