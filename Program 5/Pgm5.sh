#!/bin/bash
#SBATCH -N 32
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 00:30:00
mpicc -O3 -lm -fopenmp -mcmodel=large Pgm5.c
mpirun -np 32 ./a.out