#!/bin/bash

#./cluster_viewer start
#python -mwebbrowser http://localhost:8080
APP_ROOT="$(dirname "$(dirname "$(readlink "$0")")")"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

PYTHONPATH=$APP_ROOT $DIR/cmonkey2 "$@"

