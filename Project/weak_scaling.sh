#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -t 00:25:00
#SBATCH -p node
#SBATCH -n 64
#SBATCH -o weak_scaling_%j.txt

echo "Running Monte Carlo Sim"

mpirun --bind-to none -n 1 ./montecarlo 100000
mpirun --bind-to none -n 2 ./montecarlo 200000
mpirun --bind-to none -n 4 ./montecarlo 400000
mpirun --bind-to none -n 8 ./montecarlo 800000
mpirun --bind-to none -n 16 ./montecarlo 1600000
mpirun --bind-to none -n 32 ./montecarlo 3200000
mpirun --bind-to none -n 64 ./montecarlo 6400000



echo "All done!"