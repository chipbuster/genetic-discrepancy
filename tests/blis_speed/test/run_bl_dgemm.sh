#!/bin/bash

export DYLD_LIBRARY_PATH=/opt/intel/lib:/opt/intel/mkl/lib

#Single Thread
export OMP_NUM_THREADS=1     #Set OMP number of threads for parallel version
export BLISLAB_IC_NT=1       #Set BLISLAB number of threads for parallel version
k_start=16
k_end=1024
k_blocksize=16
echo "run_step3_st=["
echo -e "%m\t%n\t%k\t%MY_GFLOPS\t%REF_GFLOPS"
for (( k=k_start; k<=k_end; k+=k_blocksize ))
do
    ./test_bl_dgemm.x     $k $k $k 
done
echo "];"

NT=${NT:-8}

#Multi Thread
export OMP_NUM_THREADS=$NT     #Set OMP number of threads for parallel version
export BLISLAB_IC_NT=$NT       #Set BLISLAB number of threads for parallel version
k_start=64
k_end=4096
k_blocksize=64
echo "Running with $NT threads"
echo "run_step3_mt=["
echo -e "%m\t%n\t%k\t%MY_GFLOPS\t%REF_GFLOPS"
for (( k=k_start; k<=k_end; k+=k_blocksize ))
do
    ./test_bl_dgemm.x     $k $k $k 
done
echo "];"

