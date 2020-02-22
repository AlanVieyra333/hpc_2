//
//  prod_punto_vectores_mpi.c
//  
//
//  Created by Alan Vieyra on 22/02/20.
//

#include <mpi.h>
#include <stdio.h>

#define N 100

double a[N], b[N], *c;

int main(int argc, char** argv) {
    int size, rank;
    int i, l;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    c = (double*) malloc(size * sizeof(double));

    if (rank==0) {
        for (i=0; i<N; i++)
            a[i] = b[i] = i;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    l = N/size;
    
    MPI_Scatter(a, l, MPI_DOUBLE, a, l, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, l, MPI_DOUBLE, b, l, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    c[0] = 0;
    for (i=0; i<l; i++)
        c[0] += a[i]+b[i];
    
    MPI_Gather(c, 1, MPI_DOUBLE, c, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (rank==0) {
        double result = 0;
        for (i=0; i<size; i++){
            result += c[i];
            //printf("%lf\n", c[i]);
        }
        printf("Producto punto de A*B: %lf\n", result);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
