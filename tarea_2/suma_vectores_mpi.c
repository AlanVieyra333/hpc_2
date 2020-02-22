//
//  suma_vectores_mpi.c
//  
//
//  Created by Amilcar Meneses Viveros on 19/02/20.
//
//

#include <mpi.h>
#include <stdio.h>

#define N 100

double a[N], b[N], c[N];

int main(int argc, char** argv) {
    int size, rank;
    int i, l;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank==0) {
        for (i=0; i<N; i++)
            a[i] = b[i] = i;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    l = N/size;
    
    MPI_Scatter(a, l, MPI_DOUBLE, a, l, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, l, MPI_DOUBLE, b, l, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (i=0; i<l; i++)
        c[i] = a[i]+b[i];
    
    MPI_Gather(c, l, MPI_DOUBLE, c, l, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (rank==0) {
        for (i=0; i<N; i++)
            printf("%lf ", c[i]);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
