#!/bin/bash
#SBATCH -N 5
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 00:15:00
echo "Compile"
mpicc -O3 -lm Pgm3.c
echo "Execute"
mpirun -np 8 ./a.out 2147483647
echo "Done"
