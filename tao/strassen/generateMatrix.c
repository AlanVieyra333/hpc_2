/**
 * Programa para generar una matriz cuadrada de nxn.
 * Los numeros que componen la matriz estan en el rango [-100.0, 100.0].
 * Ej. Generar una matriz de 1000 x 1000:
 * ./generateMatrix.o 1000
 *
 * @author Alan Fernando Rincón Vieyra
 * @date 6/April/2020
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NUM_MAX 200.0

int main(int argc, char* argv[]) {
  if (argc != 2) {
    printf("generateMatrix <MATRIX_SIZE>\n");
    exit(1);
  }

  srand(123);

  unsigned int n = atoi(argv[1]);

  char filename[50];
  strcpy(filename, "matrix_");
  strcat(filename, argv[1]);
  strcat(filename, ".txt");

  FILE* file = fopen(filename, "w");

  for (unsigned int i = 0; i < n; i++) {
    for (unsigned int j = 0; j < n; j++) {
      if (j != 0) fprintf(file, "\t");

      fprintf(file, "%.6f", (rand() / (RAND_MAX / NUM_MAX) - (NUM_MAX / 2.0)));
    }
    fprintf(file, "\n");
  }

  fclose(file);

  return 0;
}