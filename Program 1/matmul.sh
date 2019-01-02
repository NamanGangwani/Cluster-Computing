#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node 28
#SBATCH -t 15:00:00
#SBATCH --mail-type=BEGIN
#SBATCH --mail-user=foo@bar.com
gcc -fopenmp matmul.c
./a.out