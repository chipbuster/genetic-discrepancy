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

TEST_PROGRAM="$($READLINK --canonicalize-existing "$1")"
DATADIR="$($READLINK --canonicalize-existing "$2")"
TESTSCRDIR="$(dirname "$($READLINK --canonicalize-existing "$0")")" #directory of this script

M=50 #For relatively fast tests!

POP_DIR="$(mktemp -d)"
OUT_DIR="$(mktemp -d)"

for data_file in "$DATADIR"/*; do
    echo "Running program on $data_file"
    echo  "RANDOM_SEED=0 $TEST_PROGRAM $data_file ${POP_DIR}/$(basename "$data_file").out $M"
    RANDOM_SEED=0 "$TEST_PROGRAM" "$data_file" "${POP_DIR}/$(basename "$data_file").out" "$M" &> "${OUT_DIR}/$(basename "$data_file").stdout"
done

"${TESTSCRDIR}/"test_population.sh "$POP_DIR"

if [ $? -eq 0 ]; then
    echo "Passed Population Tests"
else
    exit 1
fi

"${TESTSCRDIR}/"test_stdout.sh "$OUT_DIR"

if [ $? -eq 0 ]; then
    echo "Passed Discrepancy Tests"
    exit 0
else
    exit 1
fi
