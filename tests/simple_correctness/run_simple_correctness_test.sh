#!/usr/bin/env bash

## Dumb script to provide automated testing from CMake using test_correctness.py
if [ "$1" = "" ]; then
    echo "WARNING: this script is not meant to be run by hand!"
    exit 1
fi

TEST_PROGRAM="$(readlink --canonicalize-existing "$1")"
DATADIR="$(readlink --canonicalize-existing "$2")"
OUTPUT_DIR=/tmp/test_output
TESTSCRDIR="$(dirname "$(readlink --canonicalize-existing "$0")")" #directory of this script

M=50 #For relatively fast tests!

OUTPUT_DIR="$(mktemp -d)"

for data_file in "$DATADIR"/*; do
    echo "Running program on $data_file"
    echo  "RANDOM_SEED=0 $TEST_PROGRAM $data_file ${OUTPUT_DIR}/$(basename "$data_file").out $M"
    RANDOM_SEED=0 "$TEST_PROGRAM" "$data_file" "${OUTPUT_DIR}/$(basename "$data_file").out" "$M" &> /dev/null
done

"${TESTSCRDIR}/"test_correctness.sh "$OUTPUT_DIR"

if [ $? -eq 0 ]; then
    exit 0
else
    exit 1
fi
