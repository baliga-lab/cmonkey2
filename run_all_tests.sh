#!/bin/bash

echo "Running unit tests..."
PYTHONPATH=`pwd`/cmonkey python test/all_tests.py $@
echo "Running integration tests..."
PYTHONPATH=`pwd`/cmonkey python test/all_integration_tests.py $@
echo "Running iteration test..."
PYTHONPATH=`pwd`/cmonkey python test/iteration_test.py $@
echo "Running postprocessing test..."
PYTHONPATH=`pwd`/cmonkey python test/postproc_test.py $@
