#!/bin/bash
#SBATCH -J test-size-correct         # Job name
#SBATCH -o sizetest_output-%j.txt    # STDOUT file name
#SBATCH -e sizetest_err-%j.txt       # STDERR file name
#SBATCH -p normal                     # Partition (queue)
#SBATCH -t 03:00:00                   # Max time
#SBATCH -n 1                          # Number of MPI tasks
#SBATCH -N 1                          # Number of nodes
#SBATCH -A A-ti3                      # Allocation name
#SBATCH --mail-user=kcsong@utexas.edu
#SBATCH --mail-type=all               # Send email at begin and end of job

THREADCT=68
PROGNAME="/work/02711/chipbus/star-discrepancy-project/install-folders/projbuild/bin/test-mmult.x"
SIZES=(100 1000 10000)

NUMTRIALS=10000

export KMP_AFFINITY=balanced

for SZ in ${SIZES[@]}; do
  echo "=== Testing SIZE $SZ ==="
  "$PROGNAME" $NUMTRIALS $SIZES 300
done

