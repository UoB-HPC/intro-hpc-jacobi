//
// TODO: General description here
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 16

int main(int argc, char *argv[])
{
  // TODO: Command-line arguments for matrix size and number of iterations

  double *A    = malloc(N*N*sizeof(double));
  double *b    = malloc(N*sizeof(double));
  double *x    = malloc(N*sizeof(double));
  double *xtmp = malloc(N*sizeof(double));

  // Initialize data
  srand(0);
  for (int row = 0; row < N; row++)
  {
    double rowsum = 0.0;
    for (int col = 0; col < N; col++)
    {
      double value = rand()/(double)RAND_MAX;
      A[row + col*N] = value;
      rowsum += value;
    }
    A[row + row*N] += rowsum;
    b[row] = rand()/(double)RAND_MAX;
    x[row] = 0.0;
  }

  // Run Jacobi solver
  for (int itr = 0; itr < 10; itr++)
  {
    for (int row = 0; row < N; row++)
    {
      double tmp = 0.0;
      for (int col = 0; col < N; col++)
      {
        if (row != col)
          tmp += A[row + col*N] * x[col];
      }
      xtmp[row] = (b[row] - tmp) / A[row + row*N];
    }

    // TODO: Convergence check

    // Swap pointers
    double *tmpptr = x;
    x = xtmp;
    xtmp = tmpptr;
  }

  // Check error of final solution
  double err = 0.0;
  for (int row = 0; row < N; row++)
  {
    double tmp = 0.0;
    for (int col = 0; col < N; col++)
    {
      tmp += A[row + col*N] * x[col];
    }
    tmp = b[row] - tmp;
    err += tmp*tmp;
  }
  err = sqrt(err);
  printf("Final error = %lf\n", err);
  // TODO: Print timing information and number of iterations

  free(A);
  free(b);
  free(x);
  free(xtmp);

  return 0;
}
