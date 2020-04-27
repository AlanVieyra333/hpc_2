//
//  pi_mpi_omp.c
//
//
//  Created by Alan Vieyra on 26/04/20.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h> /* clock_t, clock, CLOCKS_PER_SEC */
#include <omp.h>
#include "mpi.h"

#define TOTAL_EXP 10000000000L

int custom_rand_next(unsigned long int *seed)  // RAND_MAX assumed to be 32767
{
  *seed = *seed * 1103515245 + 12345;
  return (unsigned int)(*seed / 65536) % RAND_MAX;
}

long getMicrotime() {
  struct timeval currentTime;
  gettimeofday(&currentTime, NULL);
  return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}

long calcula_buenos_pi(long n, int rank) {
  long i;
  double x, y;
  long success_exp = 0;

  #pragma omp parallel shared(success_exp) private(x, y)
  {
    printf("El proceso %d,%d se inicia.\n", rank, omp_get_thread_num());
    unsigned long int seed =
        getMicrotime() ^ (rank * omp_get_num_threads() + omp_get_thread_num());

    #pragma omp for
    for (i = 0; i < n; ++i) {
      x = (double)custom_rand_next(&seed) / (double)RAND_MAX;
      y = (double)custom_rand_next(&seed) / (double)RAND_MAX;
      if ((x * x + y * y) <= 1.0) {
        #pragma omp atomic
        success_exp++;
      }
    }
  }

  return success_exp;
}

int main(int argc, char **argv) {
  int size, rank;
  unsigned long eb, teb;
  int rc;
  double mi_pi;
  clock_t t1, t2;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("La tarea %d se inicia.\n", rank);

  if (rank == 0) t1 = clock();
  eb = calcula_buenos_pi(TOTAL_EXP, rank);
  rc = MPI_Reduce(&eb, &teb, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0) t2 = clock();

  if (rank == 0) {
    mi_pi = 4.0 * ((double)teb / ((double)size * TOTAL_EXP));
    printf("Aproximado de PI = %lf\n", mi_pi);
    printf("El de math.h es %lf\n", M_PI);

    printf("Time: %lf seg.\n", ((((float)t2 - (float)t1) / CLOCKS_PER_SEC)));
  }

  MPI_Finalize();
  return 0;
}
