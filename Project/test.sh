#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -t 5:00
#SBATCH -p node
#SBATCH -n 80

module load gcc openmpi



echo "Running Monte Carlo Sim"

mpirun --bind-to none -n 80 ./montecarlo 4000000




echo "All done!"