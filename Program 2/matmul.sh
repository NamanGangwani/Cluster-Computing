#!/bin/bash
#SBATCH -J matmul          # Job name
#SBATCH -o matmul.o%j      # Name of stdout output file
#SBATCH -e matmul.e%j      # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for OpenMP)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=nkg160030@utdallas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

gcc -fopenmp matmul.c

# Other commands must follow all #SBATCH directives...

module list
pwd
date

# Set thread count (default value is 1)...

export OMP_NUM_THREADS=34

# Launch OpenMP code...

./a.out        # Do not use ibrun or any other MPI launcher