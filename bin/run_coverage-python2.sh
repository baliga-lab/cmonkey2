#!/bin/bash

PYTHONPATH=`pwd` coverage2 run test/quick_tests.py $@

