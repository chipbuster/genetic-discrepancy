#!/bin/bash

SRCFILE="../../src/base/star-discrepancy.cpp"
INCDIR="../../src/include"
icc -O3 -g -p "$SRCFILE" -o knl-discrep.x -xMIC-AVX512 -std=c++11 -I "$INCDIR" -qopenmp
