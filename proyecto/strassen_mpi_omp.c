/**
 * Ej. 1000x1000
 * Time: 0.383334 seg.
 * Time with MPI: 0.820 seg.
 *
 * Edited by: Alan Fernando Rincon Vieyra
 * @date: 6/April/2020
 */

#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

// Set LEAF_SIZE to 1 if you want to the pure strassen algorithm
// otherwise, the ikj-algorithm will be applied when the split
// matrices are as small as LEAF_SIZE x LEAF_SIZE
#define LEAF_SIZE 128
#define MAXCHAR 1024 * 1024

/*
 * Implementation of the strassen algorithm, similar to
 * http://en.wikipedia.org/w/index.php?title=Strassen_algorithm&oldid=498910018#Source_code_of_the_Strassen_algorithm_in_C_language
 */

void strassen(double *A, double *B, double *C, uint32_t tam);
void divide_m(double *M, double *m11, double *m12, double *m21, double *m22,
              uint32_t tam);
void group_m(double *M, double *m11, double *m12, double *m21, double *m22,
             uint32_t tam);
void sum(double *A, double *B, double *C, uint32_t tam);
void subtract(double *A, double *B, double *C, uint32_t tam);
void readMatrix(char *filename, double *matrix);
void printMatrix(double *matrix, uint32_t tam);
void printMatrixf(double *matrix, uint32_t tam);
double *createMatrix(uint32_t tam);

uint32_t nextPowerOfTwo(uint32_t n) { return pow(2, (int)ceil(log2(n))); }

void ikjalgorithm(double *A, double *B, double *C, uint32_t tam) {
  int i, j, k;

#pragma omp parallel for private(k, j)
  for (i = 0; i < tam; i++) {
    for (k = 0; k < tam; k++) {
      for (j = 0; j < tam; j++) {
        C[i * tam + j] += A[i * tam + k] * B[k * tam + j];
      }
    }
  }
}

void strassen(double *A, double *B, double *C, uint32_t tam) {
  if (tam <= LEAF_SIZE) {
    ikjalgorithm(A, B, C, tam);
    return;
  } else {  // other cases are treated here:
    int newTam = tam / 2;
    double *a11 = createMatrix(newTam);
    double *a12 = createMatrix(newTam);
    double *a21 = createMatrix(newTam);
    double *a22 = createMatrix(newTam);
    double *b11 = createMatrix(newTam);
    double *b12 = createMatrix(newTam);
    double *b21 = createMatrix(newTam);
    double *b22 = createMatrix(newTam);
    double *c11 = createMatrix(newTam);
    double *c12 = createMatrix(newTam);
    double *c21 = createMatrix(newTam);
    double *c22 = createMatrix(newTam);
    double *p1 = createMatrix(newTam);
    double *p2 = createMatrix(newTam);
    double *p3 = createMatrix(newTam);
    double *p4 = createMatrix(newTam);
    double *p5 = createMatrix(newTam);
    double *p6 = createMatrix(newTam);
    double *p7 = createMatrix(newTam);

// dividing the matrices in 4 sub-matrices:
#pragma omp task
    { divide_m(A, a11, a12, a21, a22, tam); }
#pragma omp task
    { divide_m(B, b11, b12, b21, b22, tam); }

#pragma omp taskwait

    // Calculating p1 to p7:

#pragma omp task
    {
      double *result1 = createMatrix(newTam);
      double *result2 = createMatrix(newTam);

#pragma omp task
      {
        sum(a11, a22, result1, newTam);  // a11 + a22
      }
#pragma omp task
      {
        sum(b11, b22, result2, newTam);  // b11 + b22
      }
#pragma omp taskwait
      strassen(result1, result2, p1, newTam);  // p1 = (a11+a22) * (b11+b22)

      free(result1);
      free(result2);
    }
#pragma omp task
    {
      double *result1 = createMatrix(newTam);

      sum(a21, a22, result1, newTam);      // a21 + a22
      strassen(result1, b11, p2, newTam);  // p2 = (a21+a22) * (b11)

      free(result1);
    }
#pragma omp task
    {
      double *result1 = createMatrix(newTam);

      subtract(b12, b22, result1, newTam);  // b12 - b22
      strassen(a11, result1, p3, newTam);   // p3 = (a11) * (b12 - b22)

      free(result1);
    }
#pragma omp task
    {
      double *result1 = createMatrix(newTam);

      subtract(b21, b11, result1, newTam);  // b21 - b11
      strassen(a22, result1, p4, newTam);   // p4 = (a22) * (b21 - b11)

      free(result1);
    }
#pragma omp task
    {
      double *result1 = createMatrix(newTam);

      sum(a11, a12, result1, newTam);      // a11 + a12
      strassen(result1, b22, p5, newTam);  // p5 = (a11+a12) * (b22)

      free(result1);
    }
#pragma omp task
    {
      double *result1 = createMatrix(newTam);
      double *result2 = createMatrix(newTam);

#pragma omp task
      {
        subtract(a21, a11, result1, newTam);  // a21 - a11
      }
#pragma omp task
      {
        sum(b11, b12, result2, newTam);  // b11 + b12
      }
#pragma omp taskwait
      strassen(result1, result2, p6, newTam);  // p6 = (a21-a11) * (b11+b12)

      free(result1);
      free(result2);
    }
#pragma omp task
    {
      double *result1 = createMatrix(newTam);
      double *result2 = createMatrix(newTam);

#pragma omp task
      {
        subtract(a12, a22, result1, newTam);  // a12 - a22
      }
#pragma omp task
      {
        sum(b21, b22, result2, newTam);  // b21 + b22
      }
#pragma omp taskwait
      strassen(result1, result2, p7, newTam);  // p7 = (a12-a22) * (b21+b22)

      free(result1);
      free(result2);
    }

#pragma omp taskwait

    // calculating c21, c21, c11 e c22:
#pragma omp task
    {
      sum(p3, p5, c12, newTam);  // c12 = p3 + p5
    }
#pragma omp task
    {
      sum(p2, p4, c21, newTam);  // c21 = p2 + p4
    }
#pragma omp task
    {
      double *result1 = createMatrix(newTam);
      double *result2 = createMatrix(newTam);

      sum(p1, p4, result1, newTam);        // p1 + p4
      sum(result1, p7, result2, newTam);   // p1 + p4 + p7
      subtract(result2, p5, c11, newTam);  // c11 = p1 + p4 - p5 + p7

      free(result1);
      free(result2);
    }
#pragma omp task
    {
      double *result1 = createMatrix(newTam);
      double *result2 = createMatrix(newTam);

      sum(p1, p3, result1, newTam);        // p1 + p3
      sum(result1, p6, result2, newTam);   // p1 + p3 + p6
      subtract(result2, p2, c22, newTam);  // c22 = p1 + p3 - p2 + p6

      free(result1);
      free(result2);
    }

#pragma omp taskwait

    // Grouping the results obtained in a single matrix:
    group_m(C, c11, c12, c21, c22, tam);

    free(a11);
    free(a12);
    free(a21);
    free(a22);
    free(b11);
    free(b12);
    free(b21);
    free(b22);
    free(c11);
    free(c12);
    free(c21);
    free(c22);
    free(p1);
    free(p2);
    free(p3);
    free(p4);
    free(p5);
    free(p6);
    free(p7);
  }
}

// dividing the matrices in 4 sub-matrices.
void divide_m(double *M, double *m11, double *m12, double *m21, double *m22,
              uint32_t tam) {
  int32_t i, j;
  uint32_t newTam = tam / 2;

#pragma omp parallel for private(j)
  for (i = 0; i < newTam; i++) {
    for (j = 0; j < newTam; j++) {
      m11[i * newTam + j] = M[i * tam + j];
      m12[i * newTam + j] = M[i * tam + (j + newTam)];
      m21[i * newTam + j] = M[(i + newTam) * tam + j];
      m22[i * newTam + j] = M[(i + newTam) * tam + (j + newTam)];
    }
  }
}

// Grouping the results obtained in a single matrix.
void group_m(double *M, double *m11, double *m12, double *m21, double *m22,
             uint32_t tam) {
  int32_t i, j;
  uint32_t newTam = tam / 2;

#pragma omp parallel for private(j)
  for (i = 0; i < newTam; i++) {
    for (j = 0; j < newTam; j++) {
      M[i * tam + j] = m11[i * newTam + j];
      M[i * tam + (j + newTam)] = m12[i * newTam + j];
      M[(i + newTam) * tam + j] = m21[i * newTam + j];
      M[(i + newTam) * tam + (j + newTam)] = m22[i * newTam + j];
    }
  }
}

double *matrix_to_mPower(double *M, uint32_t tam) {
  uint32_t m = nextPowerOfTwo(tam);
  double *MPrep = createMatrix(m);
  int32_t i;

#pragma omp parallel for
  for (i = 0; i < tam; i++) {
    memcpy(&MPrep[i * m], &M[i * tam], sizeof(double) * tam);
  }

  free(M);
  return MPrep;
}

double *mPower_to_matrix(double *MPrep, uint32_t tam) {
  uint32_t m = nextPowerOfTwo(tam);
  double *M = createMatrix(tam);
  int32_t i;

#pragma omp parallel for
  for (i = 0; i < tam; i++) {
    memcpy(&M[i * tam], &MPrep[i * m], sizeof(double) * tam);
  }

  free(MPrep);
  return M;
}

void sum(double *A, double *B, double *C, uint32_t tam) {
  int32_t i, j;

#pragma omp parallel for private(j)
  for (i = 0; i < tam; i++) {
    for (j = 0; j < tam; j++) {
      C[i * tam + j] = A[i * tam + j] + B[i * tam + j];
    }
  }
}

void subtract(double *A, double *B, double *C, uint32_t tam) {
  int32_t i, j;

#pragma omp parallel for private(j)
  for (i = 0; i < tam; i++) {
    for (j = 0; j < tam; j++) {
      C[i * tam + j] = A[i * tam + j] - B[i * tam + j];
    }
  }
}

int getMatrixSize(char *filename) {
  char line[MAXCHAR];
  FILE *infile = fopen(filename, "r");
  uint32_t count = 0;

  if (infile == NULL) {
    fprintf(stderr, "Could not read file %s\n", filename);
    exit(1);
  }

  fgets(line, MAXCHAR, infile);
  fclose(infile);

  int32_t i;

  for (i = 0; line[i] != '\0'; i++) {
    if (line[i] == '\t') {
      count++;
    }
  }

  return count + 1;
}

void readMatrix(char *filename, double *matrix) {
  FILE *matrixfile = fopen(filename, "r");

  if (matrixfile == NULL) {
    fprintf(stderr, "Could not read file %s\n", filename);
    exit(1);
  }

  double num;
  int32_t i;

  for (i = 0; fscanf(matrixfile, "%lf", &num) != EOF; i++) {
    matrix[i] = num;
  }

  fclose(matrixfile);
}

void printMatrix(double *matrix, uint32_t tam) {
  uint32_t i, j;
  FILE *file = fopen("matrix_result.txt", "w");

  for (i = 0; i < tam; i++) {
    for (j = 0; j < tam; j++) {
      if (j != 0) {
        printf("\t");
      }
      printf("%.6f", matrix[i * tam + j]);
    }
    printf("\n");
  }

  fclose(file);
}

void printMatrixf(double *matrix, uint32_t tam) {
  uint32_t i, j;
  FILE *file = fopen("matrix_result.txt", "w");

  for (i = 0; i < tam; i++) {
    for (j = 0; j < tam; j++) {
      if (j != 0) {
        fprintf(file, "\t");
      }
      fprintf(file, "%.6f", matrix[i * tam + j]);
    }
    fprintf(file, "\n");
  }

  fclose(file);
}

double *createMatrix(uint32_t tam) {
  double *matrix = (double *)malloc(sizeof(double) * tam * tam);
  memset(matrix, '\0', tam * tam);
  return matrix;
}

void strassen_mpi(double *A, double *B, double *C, int n, int rank, int np) {
#pragma omp parallel
  {
    // printf("El proceso %d,%d se inicia.\n", rank, omp_get_thread_num());

#pragma omp single
    {
      strassen(A, B, C, n);  // TODO: rank
    }
  }
}

int main(int argc, char *argv[]) {
  int n;
  double *A, *B, *C;
  int tam;
  char *filename_a, *filename_b;
  int np, rank;
  double t1, t2;

  if (argc < 3) {
    printf("strassen filename_a filename_b\n");
    exit(1);
  }

  filename_a = argv[1];
  filename_b = argv[2];

  // #################### MPI INI ####################
  omp_set_nested(1);  // Enable recursive
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("La tarea %d se inicia.\n", rank);

  // Read two matrix from file only one process.
  if (rank == 0) {
    tam = getMatrixSize(filename_a);
    n = nextPowerOfTwo(tam);

    if (tam != getMatrixSize(filename_b)) {
      printf("Las matrices deben tener el mismo tamano.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    printf("Tamano: %dx%d\n", tam, tam);

    A = createMatrix(tam);
    B = createMatrix(tam);
    C = createMatrix(n);

    readMatrix(filename_a, A);
    readMatrix(filename_b, B);

    A = matrix_to_mPower(A, tam);
    B = matrix_to_mPower(B, tam);
  }

  MPI_Barrier(MPI_COMM_WORLD);  // Wait all process.

  t1 = MPI_Wtime();
  strassen_mpi(A, B, C, n, rank, np);
  t2 = MPI_Wtime();

  // Print matrix result.
  if (rank == 0) {
    printf("Time: %3.3lf seg.\n", t2 - t1);

    // Convert matrix powr of two to real size.
    C = mPower_to_matrix(C, tam);
    printMatrixf(C, tam);
  }

  free(A);
  free(B);
  free(C);

  MPI_Finalize();
  // #################### MPI END ####################

  return 0;
}
