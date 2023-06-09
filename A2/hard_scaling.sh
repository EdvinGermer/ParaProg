#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH --reservation uppmax2023-2-13_1
#SBATCH -t 5:00

echo "1 process"
mpirun --bind-to none -n 1 ./matmul ./input4.txt output.txt
echo