//
// TODO: General description here
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

static int N;
static int MAX_ITERATIONS;
static double CONVERGENCE_THRESHOLD;

// Return the current time in seconds since the Epoch
double get_timestamp();

// Parse command line arguments to set solver parameters
void parse_arguments(int argc, char *argv[]);

// Run the Jacobi solver
// Returns the number of iterations performed
int run(double *A, double *b, double *x, double *xtmp)
{
  int row, col;
  double dot;
  double diff;
  double sqdiff;
  double *ptrtmp;

  for (int itr = 0; itr < MAX_ITERATIONS; itr++)
  {
    // Perfom Jacobi iteration
    for (row = 0; row < N; row++)
    {
      dot = 0.0;
      for (col = 0; col < N; col++)
      {
        if (row != col)
          dot += A[row + col*N] * x[col];
      }
      xtmp[row] = (b[row] - dot) / A[row + row*N];
    }

    // Swap pointers
    ptrtmp = x;
    x      = xtmp;
    xtmp   = ptrtmp;

    // Check for convergence
    sqdiff = 0.0;
    for (row = 0; row < N; row++)
    {
      diff    = xtmp[row] - x[row];
      sqdiff += diff * diff;
    }
    if (sqrt(sqdiff) < CONVERGENCE_THRESHOLD)
    {
      return itr+1;
    }
  }

  return MAX_ITERATIONS;
}

int main(int argc, char *argv[])
{
  parse_arguments(argc, argv);

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
  double start = get_timestamp();
  int itr = run(A, b, x, xtmp);
  double end = get_timestamp();

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
  printf("Iterations performed = %d\n", itr);
  printf("Runtime = %lf seconds\n", (end-start));

  free(A);
  free(b);
  free(x);
  free(xtmp);

  return 0;
}

double get_timestamp()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}

int parse_int(const char *str)
{
  char *next;
  int value = strtoul(str, &next, 10);
  return strlen(next) ? -1 : value;
}

double parse_double(const char *str)
{
  char *next;
  double value = strtod(str, &next);
  return strlen(next) ? -1 : value;
}

void parse_arguments(int argc, char *argv[])
{
  // Set default values
  N = 1024;
  MAX_ITERATIONS = 100;
  CONVERGENCE_THRESHOLD = 0.0001;

  for (int i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "--convergence") || !strcmp(argv[i], "-c"))
    {
      if (++i >= argc || (CONVERGENCE_THRESHOLD = parse_double(argv[i])) < 0)
      {
        printf("Invalid convergence threshold\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--iterations") || !strcmp(argv[i], "-i"))
    {
      if (++i >= argc || (MAX_ITERATIONS = parse_int(argv[i])) < 0)
      {
        printf("Invalid number of iterations\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--norder") || !strcmp(argv[i], "-n"))
    {
      if (++i >= argc || (N = parse_int(argv[i])) < 0)
      {
        printf("Invalid matrix order\n");
        exit(1);
      }
    }
    else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
    {
      printf("\n");
      printf("Usage: ./jacobi [OPTIONS]\n\n");
      printf("Options:\n");
      printf("  -h  --help               Print this message\n");
      printf("  -c  --convergence  C     Set convergence threshold\n");
      printf("  -i  --iterations   I     Set maximum number of iterations\n");
      printf("  -n  --norder       N     Set maxtrix order\n");
      printf("\n");
      exit(0);
    }
    else
    {
      printf("Unrecognized argument '%s' (try '--help')\n", argv[i]);
      exit(1);
    }
  }
}
