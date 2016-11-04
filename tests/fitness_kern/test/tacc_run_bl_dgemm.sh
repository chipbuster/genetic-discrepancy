#!/bin/bash
#SBATCH -J bl_dgemm_job
#SBATCH -o bl_dgemm_output-%j.txt
#SBATCH -e bl_dgemm_err-%j.txt
#SBATCH -p normal
#SBATCH -t 02:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -A A-ti3
export NT=272
export KMP_AFFINITY=compact,verbose

ibrun tacc_affinity run_bl_dgemm.sh
