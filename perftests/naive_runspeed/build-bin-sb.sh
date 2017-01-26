#!/bin/bash

echo "Remember that you need to load the cxx11 module to run this."

module load cxx11/4.9.1
SRCFILE="../../src/base/star-discrepancy.cpp"
INCDIR="../../src/include"
icc -O3 -g -p "$SRCFILE" -o sb-discrep.x -xAVX -std=c++11 -I "$INCDIR" -qopenmp
