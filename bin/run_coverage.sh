#!/bin/bash

PYTHONPATH=`pwd` coverage3 run test/quick_tests.py $@

