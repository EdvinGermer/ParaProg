#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH --reservation uppmax2023-2-13_2
#SBATCH -p core -n 2
#SBATCH -t 05:00

module load gcc openmpi
mpirun --bind-to none -n 2 ./exercise1
