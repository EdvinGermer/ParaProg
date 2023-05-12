#!/bin/bash

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH -p core -n 16
#SBATCH -t 1:00:00
#SBATCH -J quicksort
#SBATCH -o strong_scaling_backwards_v2_%j.txt

echo "Running quicksort"
echo "running with pivot method 1"
mpirun --bind-to none -np 1 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 1
mpirun --bind-to none -np 2 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 1
mpirun --bind-to none -np 4 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 1
mpirun --bind-to none -np 8 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 1
mpirun --bind-to none -np 16 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 1

echo "running with pivot method 2"
mpirun --bind-to none -np 1 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 2
mpirun --bind-to none -np 2 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 2
mpirun --bind-to none -np 4 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 2
mpirun --bind-to none -np 8 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 2
mpirun --bind-to none -np 16 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 2

echo "running with pivot method 3"
mpirun --bind-to none -np 1 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 3
mpirun --bind-to none -np 2 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 3
mpirun --bind-to none -np 4 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 3
mpirun --bind-to none -np 8 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 3
mpirun --bind-to none -np 16 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/backwards125000000.txt output.txt 3

echo "All done!"