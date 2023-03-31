#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2022-2-11
#SBATCH --reservation uppmax2022-2-11_4
#SBATCH -p core -n 2
#SBATCH -t 5:00

module load gcc openmpi
mpirun -np 2 alltoall
