#!/usr/bin/env bash

## Dumb script to provide automated testing from CMake using test_population.py
if [ "$1" = "" ]; then
    echo "WARNING: this script is not meant to be run by hand!"
    exit 1
fi

if hash greadlink 2>/dev/null; then
  READLINK=greadlink
else
  READLINK=readlink
fi

TESTDIR="$1"
TESTSCRDIR="$(dirname "$($READLINK --canonicalize-existing "$0")")"
REFDIR="$TESTSCRDIR"/reference_population

"$TESTSCRDIR/test_population.py" "$TESTDIR" "$REFDIR"

if [ $? -eq 0 ]; then
    exit 0
else
    exit 1
fi
