/***************************************************************************
* File: Pgm3.c
* Author: Michael J. Quinn (original code/algorithm)
* Modification History:
*       Modified By: Naman Gangwani
*       Date Modified: 11 March  2018
* Procedures:
* main - handles and split the work amongst multiple processes to calculate
* the number of primes between 0 and the inputted number via command line
***************************************************************************/

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>

#define MIN(a, b)   ((a)<(b)?(a):(b))
#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW((id) + 1, p, n) - BLOCK_LOW(id, p, n))

/***************************************************************************
* int main( int argc, char *argv[] )
* Author: Michael J. Quinn (original code/algorithm)
* Date: 2003
* Modification History:
*       Modified By: Naman Gangwani
*       Date Modified: 11 March 2018
* Description:  handles and split the work amongst multiple processes to
* calculate the number of primes between 0 and the inputted number via
* command line. Finally, calculates sequential and parallel execution
* times as well as Karp-Flatt Metric value
*
* Parameters:
* argc I/P int The number of arguments on the command line
* argv I/P char *[] The arguments on the command line
* main O/P string Number of primes, sequential/parallel execution times,
                  Karp-Flatt Metric
* main O/P int Status code (not currently used)
**************************************************************************/
int main(int argc, char *argv[])
{
    int count;                /* Local prime count */
    double elapsed_time;      /* Parallel execution time */
    long long first;          /* Index of the first multiple */
    long long global_count;   /* Global prime count */
    long long high_value;     /* Highest value on this proc */
    int i;
    int id;                   /* Process ID number */
    long long index;          /* Index of current prime */
    long long low_value;      /* Lowest value on this proc */
    char *marked;             /* Portion of 2,...,n */
    long long n;              /* Sieving from 2,...,n */
    int p;                    /* Number of processes */
    long long proc0_size;     /* Size of proc 0's subarray */
    long long prime;          /* Current prime */
    long long size;           /* Elements in 'marked' */

    MPI_Init(&argc, &argv);

    /* Start the timer for sequential execution time */
    struct timeval start, stop; // Struct to keep track of sequential execution time
    gettimeofday(&start, NULL); // Start time before matrix multiplication

    /* Start the timer for parallel execution time*/
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* Assures that it retrieves N */
    if (argc != 2)
    {
        if (!id)
          printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoll(argv[1]); // Retrieves N from the command line

    /* Figure out this process's share of the array, as
       well as the integers represented by the first and
       last array elements  */
    low_value = 2 + BLOCK_LOW(id, p, n - 1);
    high_value = 2 + BLOCK_HIGH(id, p, n - 1);
    size = BLOCK_SIZE(id, p, n - 1);

    /* Bail out if all the primes used for sieving are
       not all held by process 0
    if ((1 + proc0_size) < (int) sqrt((double) n))
    {
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }   */

    /* Allocate this process's share of the array */
    //printf("id: %d; size: %d\n", id, size);
    marked = (char *) malloc(size);
    if (marked == NULL)
    {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    for (i = 0; i < size; i++) marked[i] = 0;
    if (!id) index = 0;
    prime = 2;
    do
    {
        /* Find the starting prime number to work with */
        if (prime * prime > low_value)
            first = prime * prime - low_value;
        else
        {
            if (!(low_value%prime))
                first = 0;
            else
                first = prime - (low_value%prime);
        }

        // Mark all multiples of k between k^2 and N
        #pragma omp parallel for
        for (i = first; i < size; i+=prime)
            marked[i] = 1; // Marked as not prime

        if (!id) // If it's the first process
        {
            while (marked[++index]); // Continue until it finds an unmarked prime number
            prime = index + 2; // Go to the index of the unmarked prime number
        }
        MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast to the first process
    } while (prime*prime <= n); // Repeat until prime*prime > n

    // Count the number of primes from its allocated size
    count = 0;
    for (i = 0; i < size; i++)
        if (!(marked[i]))
            count++;
    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); // Includes global count

    /*Stop the timer for parallel execution time */
    elapsed_time+=MPI_Wtime();

    /* Print the results if it's the first process */
    if (!id)
    {
        /* Stop the timer for sequential execution time */
        gettimeofday(&stop, NULL);
        float elapsed = ((stop.tv_sec-start.tv_sec) +
    			(stop.tv_usec-start.tv_usec)/(float)1000000); // Calculates the amount of time elapsed sequentially

        float speedup = elapsed/(elapsed_time/p); // Calculate speedup
        float e = (1/speedup + 1/p)/(1 - 1/p); // Calculate Karp-Flatt Metric
        printf("There are %d prime numbers less than or equal to %d\n", global_count, n);
        printf("Total sequential execution time: %10.6f\n", elapsed);
        printf("Total parallel execution time: %10.6f\n", elapsed_time/p);
        printf("Karp-Flatt Metric: %10.6f\n", e);
    }

    MPI_Finalize();
    return 0;
}
