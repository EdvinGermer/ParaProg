#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 4
#SBATCH -t 10:00

module load gcc openmpi

echo "64 processes"
mpirun --bind-to none -n 4 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input4.txt output.txt
echo