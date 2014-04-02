#!/bin/bash

./cluster_viewer start
python -mwebbrowser http://localhost:8080
./cmonkey.py "$@"

