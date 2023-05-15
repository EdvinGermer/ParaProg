#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -t 00:15:00
#SBATCH -p node
#SBATCH -n 64
#SBATCH -o strong_scaling_%j.txt

echo "Running Monte Carlo Sim"

mpirun --bind-to none -n 1 ./montecarlo 100000
mpirun --bind-to none -n 2 ./montecarlo 100000
mpirun --bind-to none -n 4 ./montecarlo 100000
mpirun --bind-to none -n 8 ./montecarlo 100000
mpirun --bind-to none -n 16 ./montecarlo 100000
mpirun --bind-to none -n 32 ./montecarlo 100000
mpirun --bind-to none -n 64 ./montecarlo 100000



echo "All done!"