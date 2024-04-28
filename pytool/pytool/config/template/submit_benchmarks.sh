#!/bin/bash
set -e

DIRECTORY="$1"

for subdir in "$DIRECTORY"/*; do
  if [ -d "$subdir" ]; then
    if [ -f "$subdir/submit.sh" ]; then
      pushd $subdir
      echo "Running submit.sh in $subdir"
      ./submit.sh $2
      popd
    else
      echo "submit.sh not found in $subdir"
    fi
  fi
done
