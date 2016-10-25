#!/bin/bash
#SBATCH -J bl_dgemm_job_68
#SBATCH -o bl_dgemm_output-68-%j.txt
#SBATCH -e bl_dgemm_err-68-%j.txt
#SBATCH -p normal
#SBATCH -t 02:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -A A-ti3
export OMP_NUM_THREADS=68
export BLISGEMM_IC_NT=68
export KMP_AFFINITY=balanced,verbose

ibrun tacc_affinity run_bl_dgemm.sh
