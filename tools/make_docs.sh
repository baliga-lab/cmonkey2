#!/bin/bash

for f in cmonkey/*.py ; do
  filename="${f##*/}"
  #filename="${filename%.}"
  filename=${filename%.py}
  nu_filename="$filename"
  if [ "$nu_filename" != "__init__" ]; then
    module="cmonkey.$nu_filename"
    echo $module
    pydoc -w $module
  fi
  pydoc -w cmonkey
done

