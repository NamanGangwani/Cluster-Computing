/***************************************************************************
* File: Pgm4.c
* Author: Naman Gangwani
* Date: 8 April 2018
* Procedures:
* main - handles and splits the work amongst multiple processes to
* calculate the resultant vector after multiplying a matrix with a vector
***************************************************************************/

#include <stdio.h>
#include <mpi.h>
#include <sys/time.h>

#define mpitype MPI_DOUBLE
#define M 32768
#define N 16384

#define low_value(id, p) (int) (M * id / p)
#define high_value(id, p) (int) (M/p) + low_value(id, p)

typedef double dtype;

/* Matrix, vector, and resultant vector */
dtype matrix[M][N];
dtype vec[N] = {1};
dtype outputVector[M];

/***************************************************************************
* int main( int argc, char *argv[] )
* Author: Naman Gangwani
* Date: 8 April 2018
* Description:  handles and splits the work amongst multiple processes to
* calculate the resultant vector after multiplying a matrix with a vector.
* Finally, it outputs the information and statistics of the matrix
* multiplication algorithm
*
* Parameters:
* argc I/P int      The number of arguments on the command line
* argv I/P char *[] The arguments on the command line
* main O/P string   The size of the matrix, elapsed time, and flops
* main O/P int      Status code (not currently used)
**************************************************************************/
int main(int argc, char * argv[])
{
    dtype currentOutput;     /* Dot product of current row with vector */
    int id;                   /* Process ID number */
    int p;                    /* Number of processes */
    int i, j;

    /* Populates matrix and vector with values */
    #pragma omp parallel for
    for (i = 0; i < N; i++)
      matrix[i][i] = i + 1;
    #pragma omp parallel for
    for (i = 0; i < M; i++)
      vec[i] = i + 1;

    /* Initializes MPI */
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &id); // process rank
    MPI_Comm_size(MPI_COMM_WORLD, &p); // number of processes

    /* Start the timer for execution time */
    struct timeval start, stop; // Struct to keep track of execution time
    gettimeofday(&start, NULL); // Start time before matrix-vector multiplication

    /* Splits work for this current process */
    int startBlock, end;
    startBlock = low_value(id, p);
    end = high_value(id, p);
    //printf("%d %d %d\n", startBlock, end);

    /* The calculations that this process will calculate and send */
    dtype sendBuffer[end - startBlock];
    int k = 0;

    #pragma omp parallel for private(j)
    for (i = startBlock; i < end; i++) // Loops through assigned rows
    {
        currentOutput = 0;
        #pragma omp parallel for reduction(+:currentOutput)
        for (j = 0; j < N; j++) // Goes through vector
            currentOutput+=(matrix[i][j]*vec[j]); // Takes dot product with vector

        if (!id)
          outputVector[i] = currentOutput; // Assign it to the output vector if it's the first process
        else
          sendBuffer[k++] = currentOutput; // Save in buffer to send later to process 0
        //printf("%d %d %0.2f\n", id, i, currentOutput);
    }

    // Retrieve the values calculated for the work split amongst the other processes if the current process is 0
    if(!id)
        for (i = 1; i < p; i++)
          MPI_Recv(&outputVector[low_value(i, p)], end - startBlock, mpitype, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else
        MPI_Send(&sendBuffer, end - startBlock, mpitype, 0, id, MPI_COMM_WORLD); // Send calculations over to process 0


    /* Output the information and statistics of the matrix multiplication algorithm */
    if(!id)
    {
      gettimeofday(&stop, NULL); // End time after matrix-vector multiplication
      float elapsed = ((stop.tv_sec-start.tv_sec) +
			  (stop.tv_usec-start.tv_usec)/(float)1000000); // Calculates the amount of time elapsed during matrix multiplication
      float flops = N *((float)M + (float)(M-1))/elapsed;
    	printf("M = %d, N = %d, elapsed = %g, flops = %g\n", M, N, elapsed, flops);
    }

//    if (!id)
//        for (i = 0; i < M; i++)
//            printf("%0.2f\n", outputVector[i]);

    MPI_Finalize();
    return 0;
}

