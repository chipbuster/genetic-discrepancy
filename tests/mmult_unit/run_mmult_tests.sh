#!/bin/bash
#SBATCH -J test-discrep-speed         # Job name
#SBATCH -o speedtest_output-%j.txt    # STDOUT file name
#SBATCH -e speedtest_err-%j.txt       # STDERR file name
#SBATCH -p normal                     # Partition (queue)
#SBATCH -t 02:00:00                   # Max time
#SBATCH -n 1                          # Number of MPI tasks
#SBATCH -N 1                          # Number of nodes
#SBATCH -A A-ti3                      # Allocation name
#SBATCH --mail-user=kcsong@utexas.edu
#SBATCH --mail-type=all               # Send email at begin and end of job

THREADCT=68
PROGNAME=true

export KMP_AFFINITY=balanced

# Your commands here!
