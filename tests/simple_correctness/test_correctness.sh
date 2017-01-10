#!/usr/bin/env bash

## Dumb script to provide automated testing from CMake using test_correctness.py
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
REFDIR="$TESTSCRDIR"/reference_output

"$TESTSCRDIR/test_correctness.py" "$TESTDIR" "$REFDIR"

if [ $? -eq 0 ]; then
    exit 0
else
    exit 1
fi
