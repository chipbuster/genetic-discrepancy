#!/bin/bash
#SBATCH -J test-discrep-speed
#SBATCH -o speedtest_output-%j.txt
#SBATCH -e speedtest_err-%j.txt
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -A A-ti3

MVALUES=(1000 10000 100000)
THREADCT=16
PROGNAME=sb-discrep.x

module load cxx11/4.9.1

# There should really only be one datafile anyways
DATAFILE=$(find . -name 'pts*' | head -n 1)

export KMP_AFFINITY=compact

for M in ${MVALUES[@]}; do
  echo "Running $M"
  time ./$PROGNAME $DATAFILE $THREADCT ${M}.out &> ${M}-stdout.log
  echo "-------"
done
