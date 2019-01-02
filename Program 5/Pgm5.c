/***************************************************************************
* File: Pgm5.c
* Author: Naman Gangwani
* Date: 1 May 2018
* Procedures:
* main - handles and splits the work amongst multiple processes to
* sort the 32-bit numbers with communication between them, ultimately
* implementing the PSRS algorithm
* quicksort - recursively sorts a given list by splitting them
* based on selected pivot values
* swap - swaps the contents of the two variables at their memory addresses
* printArray - prints the elements of a given array in a single line
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>

#define mpitype MPI_INT
#define N 128000000

void quicksort(int *arr, int low, int high);
void swap(int* a, int* b);
void printArray(int *arr, int size);

/***************************************************************************
* int main( int argc, char *argv[] )
* Author: Naman Gangwani
* Date: 1 May 2018
* Description: handles and splits the work amongst multiple processes to
* sort the 32-bit numbers with communication between them, ultimately
* implementing the PSRS algorithm
*
* Parameters:
* argc I/P int      The number of arguments on the command line
* argv I/P char *[] The arguments on the command line
* main O/P string   Elapsed time when sorting and samples per second
* main O/P int      Status code (not currently used)
**************************************************************************/
int main(int argc, char * argv[])
{
    int id;                   /* Process ID number */
    int p;                    /* Number of processes */
    int i, j;
    int status, gstatus;
    time_t t;
   
    /* Initializes MPI */
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &id); // process rank
    MPI_Comm_size(MPI_COMM_WORLD, &p); // number of processes
    
    /* Intializes random number generator */
    srand((unsigned) time(&t) * id);
    
    /* Creates array and randomly generates numbers in it to fill it */
  	int *arr = (int *)malloc(N * sizeof(int));
  	#pragma omp parallel for
  	for (i = 0; i < N; i++)
      arr[i] = (int) rand()%500;
    
    /* Start the timer for execution time of sorting */
    struct timeval start, stop; // Struct to keep track of execution time
    gettimeofday(&start, NULL); // Start time before matrix-vector multiplication
    
    /* Quicksorts the array with randomly generated number */
    //printf("Process %d initial: ", id); fflush (stdout); printArray(arr, N);
  	quicksort(arr, 0, N - 1);
    
    /* Selectrs regular samples */
    int regularSamples[p];
    j = 0;
    for (i = N/p; i < N; i+=(N/p))
        regularSamples[j++] = arr[i];

    int pivots[p-1]; // Pivots to retrieve from calculations by one process
    if (id != 0) // If it's not process 0
    {
        /* Send regular values, retrieve pivot values after they are calculated */
        MPI_Send(&regularSamples, p, mpitype, 0, id, MPI_COMM_WORLD); // Send regular samples over to process 0
        MPI_Recv(&pivots, p-1, mpitype, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        /* Populates array with regular samples form all processes */
        int *toSort = (int *)malloc(p*p*sizeof(int));
        #pragma omp parallel for
        for (i = 0; i < p; i++)
            toSort[i] = regularSamples[i];
        // Receive regular sample from other processes 
        for (i = 1; i < p; i++)
            MPI_Recv(&toSort[i*p], p, mpitype, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Sort the regular samples
        quicksort(toSort, 0, p*p-1);
        
        // Calculate the pivots
        j = 0;
        int inc = p*p/(p-1);
        for (i = inc - inc/2; i < p*p; i+=inc)
            pivots[j++] = toSort[i];
            
        // Send over pivots to all other processes
        for (i = 1; i < p; i++)
            MPI_Send(&pivots, p-1, mpitype, i, 0, MPI_COMM_WORLD);
        free(toSort);
    }
  
    /* Calculate the send count and the send displacements to send to each process based on the pivots */
    int lastIndex = 0;
    int *sendCount = (int *)malloc(p * sizeof(int));
    int *sendDispl = (int *)malloc(p * sizeof(int));
    //#pragma omp parallel for private(j)
    for (i = 0; i < p-1; i++)
    {
        sendCount[i] = 0;
        sendDispl[i] = lastIndex;
        for (j = lastIndex; j < N; j++) // for each process
            if (arr[j] <= pivots[i]) // If it is within the pivot value, increment
                sendCount[i] = sendCount[i] + 1;
            else
            {
                // Done with current process, check next process
                lastIndex = j;
                break;
            }
    }
    sendCount[p-1] = N - lastIndex;
    sendDispl[p-1] = lastIndex;
    
    /* Implementation of MPI_Alltoall to retrieve the receiveCount and sendCount for each process */
    int *receiveCount = (int *)malloc(p * sizeof(int));
    int *receiveDispl = (int *)malloc(p * sizeof(int));
    // Send the send count and send displacement from all other processes
    for (i = 0; i < p; i++)
      if (i != id)
      {
        MPI_Send(&sendCount[i], 1, mpitype, i, id, MPI_COMM_WORLD);
        MPI_Send(&sendDispl[i], 1, mpitype, i, id, MPI_COMM_WORLD);
      }
    // Receive the receive count and receive displacement from all other processes
    for (i = 0; i < p; i++)
      if (i != id)
      {
        MPI_Recv(&receiveCount[i], 1, mpitype, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&receiveDispl[i], 1, mpitype, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      else
      {
        receiveCount[i] = sendCount[i];
        receiveDispl[i] = sendDispl[i];
      }
    
    /* Calculate the total amount needed to allocate in memory when receiving buffers of partitions */
    int total = 0;
    #pragma omp parallel for reduction(+:total)
    for (i = 0; i < p; i++)
        total+=receiveCount[i];
    int *newArr = malloc(total * sizeof(int));
    
    /* Send partitions based on send counts and send displacements to all other processes */
    int displ = 0;
    for (i = 0; i < p; i++)
      if (i != id)
        MPI_Send(&arr[sendDispl[i]], sendCount[i], mpitype, i, id, MPI_COMM_WORLD);
    /* Populates final array with the correct partitions to the process */
    for (i = 0; i < p; i++)
    {
      // Receives partition from each process and populates new array with it
      if (i != id)
      {
        MPI_Recv(&newArr[displ], receiveCount[i], mpitype, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        displ+=receiveCount[i];
      }
      else
        for (j = receiveDispl[i]; j < receiveDispl[i] + receiveCount[i]; j++)
          newArr[displ++] = arr[j];
    }
    
    /* Free uncessesary arrays from memory */
    free(arr);
    free(sendCount);
    free(sendDispl);
    free(receiveCount);
    free(receiveDispl);
    
    /* Sorts the process's sample to get the final result */
    quicksort(newArr, 0, total - 1);
    //printf("Process %d final: ", id); fflush (stdout); printArray(newArr, total);
    
    /* Free final array from memory*/
    free(newArr);
    
    /* Time elapsed */
    if (!id)
    {
      gettimeofday(&stop, NULL); // End time after matrix-vector multiplication
      float elapsed = ((stop.tv_sec-start.tv_sec) +
  			  (stop.tv_usec-start.tv_usec)/(float)1000000); // Calculates the amount of time elapsed during PSRS
      printf("\nElapsed time for PSRS algorithm to complete: %0.2f seconds\n", elapsed); fflush (stdout)
      float samplesPerSec = (float)(p*N)/elapsed;
      printf("Samples per second: %0.2f\n", samplesPerSec); fflush (stdout)
    }
                                                        
    
    /* Finalize and exit */
  	MPI_Finalize();
    return 0;
}

/***************************************************************************
* void quicksort( int *arr, int low, int high )
* Author: Naman Gangwani
* Date: 1 May 2018
* Description: recursively sorts a given list by splitting them
* based on selected pivot values 
*
* Parameters:
* arr  I/P int *     The array to be sorted
* low  I/P int       The low value to begin sorting from
* high I/P int       The high value to end the recursive sorting
**************************************************************************/
void quicksort(int *arr, int low, int high) {
    if (low >= high) // Cancel because sort is complete
        return;
    
    // Calculate pivot
    int pivot = low;
    int last = pivot;
    int i;
    
    // Swaps back and forth from low value up to high value
    for (i = pivot + 1; i <= high; i++)
        if (arr[i] <= arr[pivot])
          swap(&arr[++last], &arr[i]);

    // Swap last one to be higher than the pivot
    swap(&arr[last], &arr[pivot]);
    
    // Continue sorting recursively
    #pragma omp parallel sections
    {
      #pragma omp section
      {
        quicksort(arr, low, last - 1);
      }
      #pragma omp section
      {
        quicksort(arr, last + 1, high);
      }
    }
}

/***************************************************************************
* int swap( int* a, int *b )
* Author: Naman Gangwani
* Date: 1 May 2018
* Description: recursively sorts a given list by splitting them
* based on selected pivot values 
*
* Parameters:
* a  I/P int*       The first value to be swapped with the second one
* b  I/P int*       The second value to be swapped with the first one
**************************************************************************/
void swap(int* a, int* b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

/***************************************************************************
* int swap( int* a, int *b )
* Author: Naman Gangwani
* Date: 1 May 2018
* Description: prints the elements of a given array in a single line
*
* Parameters:
* arr  I/P int*       The array whose elements is to be printed
* size I/P int        The size of the array
*      O/P string     The elements of the given array prnted on a single
                      line followed by a new line character
**************************************************************************/
void printArray(int *arr, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        printf("%ld ", arr[i]);
        fflush (stdout);
    }
	  printf("\n");
    fflush (stdout);
}
