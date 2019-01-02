///////////////////////////////////////////////////////////////////////////////
// matmul.c
// Author: Dr. Richard A. Goodrum, Ph.D.
//
// Procedures:
//
// main	generates matrices and tests matmul
// matmul	basic, brute force matrix multiply
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <sys/time.h>
#include <omp.h>

int matmul(int l, int m, int n, float *A, float *B, float *C);

///////////////////////////////////////////////////////////////////////////////
// int main(int argc, char *argv[])
// Author: Dr. Richard A. Goodrum, Ph.D.
// Date:  16September 2017
// Description: Generates two matrices and then calls matmul to multiply them.
// 	Finally, it verifies that the results are correct.
//
// Parameters:
//	argc	I/P	int	The number of arguments on the command line
//	argv	I/P	char *[]	The arguments on the command line
//	main	O/P	int	Status code
///////////////////////////////////////////////////////////////////////////////
#ifndef L
#define L (250*1024/256)
#endif
#ifndef M
#define M (250*1024/256)
#endif
#ifndef N
#define N (250*1024/256)
#endif
//#define DEBUG 1
float A[L*M], B[M*N], C[L*N];

int main(int argc, char *argv[])
{
	int i, j, k;

  /* Populates A array to an upper triangular matrix of various numbers*/
	for(i = 0; i < L; i++) // Loops through A's first dimension
		for(j = 0; j < M; j++) // Loops through A's second dimension
			if(i <= j) // If it is in the upper triangle or diagonal
			{
				A[i*M + j] = (float) (i*M + j + 1); // Compute a value for it
			}
			else
			{
				A[i*M + j] = 0.0; // Otherwise, it is 0 in lower triangle
				//A[i*M+j] = (float) (i*M+j+1);
			}

  /* Populates B array iteratively */
	for(j = 0; j < M; j++) // Loops through B's first dimension
		for(k = 0; k < N; k++) // Loops through B's second dimension
			if(j <= k) // If it's in the upper triangle or diagonal
			{
			  if(k < M) // If it is a part of the square part of the matrix
				B[j*N + k] = 1.0; // Set it to 1
			  else
				//Otherwise, set it to be one plus the previous value of its row
				B[j*N + k] = B[j*N + k - 1] + 1.0;
			}
			else
			{
				B[j*N + k] = 0.0; // Otherwise, it is 0 in lower triangle
			}

  /* Populates C array  with temporary numbers*/
	for(i = 0; i < L; i++) // Loops through C's first dimension
	  for(k = 0; k < N; k++) // Loops through C's second dimension
	    C[i*N + k] = - (float) L*M*N; // Assigns it to an arbritrary value

	struct timeval start, stop; // Struct to keep track of execution time
	gettimeofday(&start, NULL); // Start time before matrix multiplication
	matmul(L, M, N, A, B, C);
	gettimeofday(&stop, NULL); // End time after matrix multiplication
	float elapsed = ((stop.tv_sec-start.tv_sec) +
			(stop.tv_usec-start.tv_usec)/(float)1000000); // Calculates the amount of time elapsed during matrix multiplication
	
	// Number of floating point operations per second that occurred during matrix multiplication
	float flops = (2 * (float)L * (float)M * (float)N)/elapsed; 
	
	// Outputs the information and statistics of the matrix multiplication algorithm
	printf("L = %d, M = %d, N = %d, elapsed = %g, flops = %g\n",
		L, M, N, elapsed, flops);
		
    printf("\n");

/* Prints matrices to check for correctness only if DEBUG mode is on */
#ifdef DEBUG
	/* Print the A matrix */
	printf( "A:\n" );
	for(i = 0; i < L; i++) // For each row
	{
	  printf("%g", A[i*M]); 
	  for(j = 1; j < M; j++) // For each column of each row
	  {
	    printf(" %g", A[i*M + j]); // Print value at that position
	  }
	  printf("\n");
	}
	
	/* Print the B matrix */
	printf("B:\n");
	for(j = 0; j < M; j++) // For each row
	{
	  printf("%g", B[j*N]);
	  for(k = 1; k < N; k++) // For each column of each row
	  {
	    printf(" %g", B[j*N + k]); // Print value at that position
	  }
	  printf("\n");
	}
	
	/* Print the C matrix */
	printf("C:\n");
	for(i = 0; i < L; i++)
	{
	  printf("%g", C[i*N]);
	  for(k = 1; k < N; k++) // For each row
	  {
	    printf(" %g", C[i*N + k]); // For each column of each row
	  }
	  printf("\n");
	}
#endif
}

///////////////////////////////////////////////////////////////////////////////
// int main(int argc, char *argv[])
// Author: Dr. Richard A. Goodrum, Ph.D.
// Date:  16September 2017
// Modified by: Naman Gangwani
// Date: 25February 2018
// Description: Multiplies two matrices 
//              and stores the results into a third matrix.
//
// Parameters:
//	l	I/P	int	The first dimension of A and C
//	m	I/P	int	The second dimension of A and  first of B
//	n	I/P	int	The second dimension of B and C
//	A	I/P	float *	The first input matrix
//	B	I/P	float *	The second input matrix
//	C	O/P	float *	The output matrix
//	matmul	O/P	int	Status code
///////////////////////////////////////////////////////////////////////////////
int matmul(int l, int m, int n, float *A, float *B, float *C)
{
	int i, j, k, result;
	
  #pragma omp parallel for private(j, k)
  for(i = 0; i < l; i++) // For each row of A and C
    for(k = 0; k < n; k++) // For each column of B and C
    {
      result = 0;
      /* Take the dot product of the row of A with the column of B */
      #pragma omp parallel for reduction(+:result)
      for(j = 0; j < m; j++) // For each column of A
      {
        // Add single component of dot product
          result += A[i*m + j] * B[j*n + k];
      }
      C[i*n + k] = result; // Save the dot product into C
    }
}
