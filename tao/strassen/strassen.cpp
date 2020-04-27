/**
 * Ej. 1000x1000
 * Time without pointers: 13 seg.
 * Time without memcpy nor -O2: 2.528321 seg.
 * Time: 0.394698 seg.
 *
 * Edited by: Alan Fernando Rincon Vieyra
 * @date: 6/April/2020
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

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
uint32_t nextPowerOfTwo(uint32_t n);
void strassenR(double *A, double *B, double *C, uint32_t tam);
void sum(double *A, double *B, double *C, uint32_t tam);
void subtract(double *A, double *B, double *C, uint32_t tam);
void read(char *filename, double *matrix);
void printMatrix(double *matrix, uint32_t tam);
void printMatrixf(double *matrix, uint32_t tam);
double *createMatrix(uint32_t tam);

uint32_t nextPowerOfTwo(uint32_t n) { return pow(2, int(ceil(log2(n)))); }

void ikjalgorithm(double *A, double *B, double *C, uint32_t tam) {
  for (int i = 0; i < tam; i++) {
    for (int k = 0; k < tam; k++) {
      for (int j = 0; j < tam; j++) {
        C[i * tam + j] += A[i * tam + k] * B[k * tam + j];
      }
    }
  }
}

void strassenR(double *A, double *B, double *C, uint32_t tam) {
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
    double *aResult = createMatrix(newTam);
    double *bResult = createMatrix(newTam);

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

    sum(a11, a22, aResult, newTam);           // a11 + a22
    sum(b11, b22, bResult, newTam);           // b11 + b22
    strassenR(aResult, bResult, p1, newTam);  // p1 = (a11+a22) * (b11+b22)

    sum(a21, a22, aResult, newTam);       // a21 + a22
    strassenR(aResult, b11, p2, newTam);  // p2 = (a21+a22) * (b11)

    subtract(b12, b22, bResult, newTam);  // b12 - b22
    strassenR(a11, bResult, p3, newTam);  // p3 = (a11) * (b12 - b22)

    subtract(b21, b11, bResult, newTam);  // b21 - b11
    strassenR(a22, bResult, p4, newTam);  // p4 = (a22) * (b21 - b11)

    sum(a11, a12, aResult, newTam);       // a11 + a12
    strassenR(aResult, b22, p5, newTam);  // p5 = (a11+a12) * (b22)

    subtract(a21, a11, aResult, newTam);      // a21 - a11
    sum(b11, b12, bResult, newTam);           // b11 + b12
    strassenR(aResult, bResult, p6, newTam);  // p6 = (a21-a11) * (b11+b12)

    subtract(a12, a22, aResult, newTam);      // a12 - a22
    sum(b21, b22, bResult, newTam);           // b21 + b22
    strassenR(aResult, bResult, p7, newTam);  // p7 = (a12-a22) * (b21+b22)

    // calculating c21, c21, c11 e c22:

    sum(p3, p5, c12, newTam);  // c12 = p3 + p5
    sum(p2, p4, c21, newTam);  // c21 = p2 + p4

    sum(p1, p4, aResult, newTam);        // p1 + p4
    sum(aResult, p7, bResult, newTam);   // p1 + p4 + p7
    subtract(bResult, p5, c11, newTam);  // c11 = p1 + p4 - p5 + p7

    sum(p1, p3, aResult, newTam);        // p1 + p3
    sum(aResult, p6, bResult, newTam);   // p1 + p3 + p6
    subtract(bResult, p2, c22, newTam);  // c22 = p1 + p3 - p2 + p6

    // Grouping the results obtained in a single matrix:
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
  }
}

void strassen(double *A, double *B, double *C, uint32_t tam) {
  // uint32_t m = tam;
  uint32_t m = nextPowerOfTwo(tam);
  double *APrep = createMatrix(m);
  double *BPrep = createMatrix(m);
  double *CPrep = createMatrix(m);

  // TODO: implement memcpy
  for (uint32_t i = 0; i < tam; i++) {
    memcpy(&APrep[i * m], &A[i * tam], sizeof(double) * tam);
    memcpy(&BPrep[i * m], &B[i * tam], sizeof(double) * tam);
  }

  strassenR(APrep, BPrep, CPrep, m);

  for (uint32_t i = 0; i < tam; i++) {
    memcpy(&C[i * tam], &CPrep[i * m], sizeof(double) * tam);
  }

  free(APrep);
  free(BPrep);
  free(CPrep);
}

void sum(double *A, double *B, double *C, uint32_t tam) {
  uint32_t i, j;

  for (i = 0; i < tam; i++) {
    for (j = 0; j < tam; j++) {
      C[i * tam + j] = A[i * tam + j] + B[i * tam + j];
    }
  }
}

void subtract(double *A, double *B, double *C, uint32_t tam) {
  uint32_t i, j;

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

  for (uint32_t i = 0; line[i] != '\0'; i++) {
    if (line[i] == '\t') {
      count++;
    }
  }

  return count + 1;
}

void read(char *filename, double *matrix, uint32_t n) {
  char line[MAXCHAR];
  FILE *matrixfile = fopen(filename, "r");

  if (matrixfile == NULL) {
    fprintf(stderr, "Could not read file %s\n", filename);
    exit(1);
  }

  double a;
  for (uint32_t i = 0; fgets(line, MAXCHAR, matrixfile) != NULL; i++) {
    std::istringstream iss(line);
    for (uint32_t j = 0; iss >> matrix[i * n + j]; j++);
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

int main(int argc, char *argv[]) {
  char *filename_a, *filename_b;

  if (argc != 3) {
    printf("strassen filename_a filename_b\n");
    exit(1);
  }

  filename_a = argv[1];
  filename_b = argv[2];

  int n = getMatrixSize(filename_a);

  if (n != getMatrixSize(filename_b)) {
    printf("Las matrices deben tener el mismo tamano.\n");
    exit(1);
  }

  printf("Tamano: %dx%d\n", n, n);

  double *A = createMatrix(n);
  double *B = createMatrix(n);
  double *C = createMatrix(n);

  read(filename_a, A, n);
  read(filename_b, B, n);

  clock_t t1, t2;
  t1 = clock();
  strassen(A, B, C, n);
  t2 = clock();

  printf("Time: %lf seg.\n", ((((float)t2 - (float)t1) / CLOCKS_PER_SEC)));

  printMatrixf(C, n);

  free(A);
  free(B);
  free(C);

  return 0;
}