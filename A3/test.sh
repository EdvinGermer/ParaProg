#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 10:00

module load gcc openmpi

echo "16 processes"
mpirun --bind-to none -n 16 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input4.txt output.txt 1
echo