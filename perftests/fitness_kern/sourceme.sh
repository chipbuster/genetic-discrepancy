#!/bin/bash
export BLISLAB_DIR=.
echo "BLISLAB_DIR = $BLISLAB_DIR"

# Compiler Options (true: Intel compiler; false: GNU compilers)
export BLISLAB_USE_INTEL=true
echo "BLISLAB_USE_INTEL = $BLISLAB_USE_INTEL"

# Whether reference implementation uses BLAS or not?
export BLISLAB_USE_BLAS=false
echo "BLISLAB_USE_BLAS = $BLISLAB_USE_BLAS"

# Optimization Level (O0, O1, O2, O3)
export COMPILER_OPT_LEVEL=O3
echo "COMPILER_OPT_LEVEL = $COMPILER_OPT_LEVEL"

# Manually set the BLAS path if BLIS_USE_BLAS=true and using GNU compiler.
export BLAS_DIR=$WORK/star-discrepancy-project/install-folders/blis-from-flame/
echo "BLAS_DIR = $BLAS_DIR"

# Parallel Options
export KMP_AFFINITY=compact,verbose
