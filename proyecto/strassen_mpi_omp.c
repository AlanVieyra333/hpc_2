/**
 * Ej. 1000x1000
 * Time: 0.383334 seg.
 *
 * Edited by: Alan Fernando Rincon Vieyra
 * @date: 26/April/2020
 */

#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

double *strassen(double *A, double *B, uint32_t tam);
void preproc_strassen(double *A, double *B, double *C, uint32_t tam);
uint32_t nextPowerOfTwo(uint32_t n);
double *sum(double *A, double *B, uint32_t tam);
double *subtract(double *A, double *B, uint32_t tam);
double *createEmpyMatrix(uint32_t tam);
void readMatrix(char *filename, double *matrix);
void printMatrix(double *matrix, uint32_t tam);
void printMatrixf(double *matrix, uint32_t tam);

uint32_t nextPowerOfTwo(uint32_t n) { return pow(2, (int)ceil(log2(n))); }

double *ikjalgorithm(double *A, double *B, uint32_t tam) {
  double *C = createEmpyMatrix(tam);

  for (int i = 0; i < tam; i++) {
    for (int k = 0; k < tam; k++) {
      for (int j = 0; j < tam; j++) {
        C[i * tam + j] += A[i * tam + k] * B[k * tam + j];
      }
    }
  }

  return C;
}

double *strassen(double *A, double *B, uint32_t tam) {
  if (tam <= LEAF_SIZE) {
    double *C = ikjalgorithm(A, B, tam);
    return C;
  } else {  // other cases are treated here:
    int newTam = tam / 2;
    double *p1, *p2, *p3, *p4, *p5, *p6, *p7;
    double *a11 = createEmpyMatrix(newTam);
    double *a12 = createEmpyMatrix(newTam);
    double *a21 = createEmpyMatrix(newTam);
    double *a22 = createEmpyMatrix(newTam);

    double *b11 = createEmpyMatrix(newTam);
    double *b12 = createEmpyMatrix(newTam);
    double *b21 = createEmpyMatrix(newTam);
    double *b22 = createEmpyMatrix(newTam);

    double *c11 = createEmpyMatrix(newTam);
    double *c12 = createEmpyMatrix(newTam);
    double *c21 = createEmpyMatrix(newTam);
    double *c22 = createEmpyMatrix(newTam);

    double *aResult = createEmpyMatrix(newTam);
    double *bResult = createEmpyMatrix(newTam);

    // dividing the matrices in 4 sub-matrices:
    uint32_t i, j;
    for (i = 0; i < newTam; i++) {
      for (j = 0; j < newTam; j++) {
        a11[i * newTam + j] = A[i * tam + j];
        a12[i * newTam + j] = A[i * tam + (j + newTam)];
        a21[i * newTam + j] = A[(i + newTam) * tam + j];
        a22[i * newTam + j] = A[(i + newTam) * tam + (j + newTam)];

        b11[i * newTam + j] = B[i * tam + j];
        b12[i * newTam + j] = B[i * tam + (j + newTam)];
        b21[i * newTam + j] = B[(i + newTam) * tam + j];
        b22[i * newTam + j] = B[(i + newTam) * tam + (j + newTam)];
      }
    }

    // Calculating p1 to p7:

    aResult = sum(a11, a22, newTam);          // a11 + a22
    bResult = sum(b11, b22, newTam);          // b11 + b22
    p1 = strassen(aResult, bResult, newTam);  // p1 = (a11+a22) * (b11+b22)

    aResult = sum(a21, a22, newTam);      // a21 + a22
    p2 = strassen(aResult, b11, newTam);  // p2 = (a21+a22) * (b11)

    bResult = subtract(b12, b22, newTam);  // b12 - b22
    p3 = strassen(a11, bResult, newTam);   // p3 = (a11) * (b12 - b22)

    bResult = subtract(b21, b11, newTam);  // b21 - b11
    p4 = strassen(a22, bResult, newTam);   // p4 = (a22) * (b21 - b11)

    aResult = sum(a11, a12, newTam);      // a11 + a12
    p5 = strassen(aResult, b22, newTam);  // p5 = (a11+a12) * (b22)

    aResult = subtract(a21, a11, newTam);     // a21 - a11
    bResult = sum(b11, b12, newTam);          // b11 + b12
    p6 = strassen(aResult, bResult, newTam);  // p6 = (a21-a11) * (b11+b12)

    aResult = subtract(a12, a22, newTam);     // a12 - a22
    bResult = sum(b21, b22, newTam);          // b21 + b22
    p7 = strassen(aResult, bResult, newTam);  // p7 = (a12-a22) * (b21+b22)

    // calculating c21, c21, c11 e c22:

    c12 = sum(p3, p5, newTam);  // c12 = p3 + p5
    c21 = sum(p2, p4, newTam);  // c21 = p2 + p4

    aResult = sum(p1, p4, newTam);        // p1 + p4
    bResult = sum(aResult, p7, newTam);   // p1 + p4 + p7
    c11 = subtract(bResult, p5, newTam);  // c11 = p1 + p4 - p5 + p7

    aResult = sum(p1, p3, newTam);        // p1 + p3
    bResult = sum(aResult, p6, newTam);   // p1 + p3 + p6
    c22 = subtract(bResult, p2, newTam);  // c22 = p1 + p3 - p2 + p6

    // Grouping the results obtained in a single matrix:
    double *C = createEmpyMatrix(tam);
    for (i = 0; i < newTam; i++) {
      for (j = 0; j < newTam; j++) {
        C[i * tam + j] = c11[i * newTam + j];
        C[i * tam + (j + newTam)] = c12[i * newTam + j];
        C[(i + newTam) * tam + j] = c21[i * newTam + j];
        C[(i + newTam) * tam + (j + newTam)] = c22[i * newTam + j];
      }
    }

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
    free(aResult);
    free(bResult);

    return C;
  }
}

double *sum(double *A, double *B, uint32_t tam) {
  double *C = createEmpyMatrix(tam);
  uint32_t i, j;

  for (i = 0; i < tam; i++) {
    for (j = 0; j < tam; j++) {
      C[i * tam + j] = A[i * tam + j] + B[i * tam + j];
    }
  }

  return C;
}

double *subtract(double *A, double *B, uint32_t tam) {
  double *C = createEmpyMatrix(tam);
  uint32_t i, j;

  for (i = 0; i < tam; i++) {
    for (j = 0; j < tam; j++) {
      C[i * tam + j] = A[i * tam + j] - B[i * tam + j];
    }
  }

  return C;
}

/**
 * Funcion que preprocesa las matrices para que estas tengan un tamano
 * igual a una potencia de 2 (la inmediata superior).
 */
void preproc_strassen(double *A, double *B, double *C, uint32_t tam) {
  // uint32_t m = tam;
  uint32_t m = nextPowerOfTwo(tam);
  double *APrep = createEmpyMatrix(m);
  double *BPrep = createEmpyMatrix(m);
  double *CPrep;

  for (uint32_t i = 0; i < tam; i++) {
    memcpy(&APrep[i * m], &A[i * tam], sizeof(double) * tam);
    memcpy(&BPrep[i * m], &B[i * tam], sizeof(double) * tam);
  }

  CPrep = strassen(APrep, BPrep, m);

  for (uint32_t i = 0; i < tam; i++) {
    memcpy(&C[i * tam], &CPrep[i * m], sizeof(double) * tam);
  }

  free(APrep);
  free(BPrep);
  free(CPrep);
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

  for (uint32_t i = 0; line[i] != '\0'; i++) {
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
  for (uint32_t i = 0; fscanf(matrixfile, "%lf", &num) != EOF; i++) {
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

double *createEmpyMatrix(uint32_t tam) {
  double *matrix = (double *)malloc(sizeof(double) * tam * tam);
  memset(matrix, '\0', tam * tam);
  return matrix;
}

int main(int argc, char *argv[]) {
  int n;
  double *A, *B, *C;
  char *filename_a, *filename_b;
  int size, rank;
  clock_t t1, t2;

  if (argc < 3) {
    printf("strassen filename_a filename_b\n");
    exit(1);
  }

  filename_a = argv[1];
  filename_b = argv[2];

  n = getMatrixSize(filename_a);

  if (n != getMatrixSize(filename_b)) {
    printf("Las matrices deben tener el mismo tamano.\n");
    exit(1);
  }

  printf("Tamano: %dx%d\n", n, n);

  A = createEmpyMatrix(n);
  B = createEmpyMatrix(n);
  C = createEmpyMatrix(n);

  readMatrix(filename_a, A);
  readMatrix(filename_b, B);

  // #################### MPI INI ####################

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("La tarea %d se inicia.\n", rank);

  if (rank == 0) t1 = clock();
  preproc_strassen(A, B, C, n);  // TODO: rank
  if (rank == 0) t2 = clock();

  if (rank == 0) {
    printf("Time: %lf seg.\n", ((((float)t2 - (float)t1) / CLOCKS_PER_SEC)));
    printMatrixf(C, n);
  }

  MPI_Finalize();

  // #################### MPI END ####################

  free(A);
  free(B);
  free(C);

  return 0;
}
