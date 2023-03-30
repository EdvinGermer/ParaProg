#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH --reservation uppmax2023-2-13_1
#SBATCH -t 5:00

module load gcc openmpi
mpirun integral2d
