#!/bin/bash

PYTHONPATH=`pwd` coverage run test/quick_tests.py $@

