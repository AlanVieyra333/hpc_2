//
//  status.c
//  
//
//  Created by Amilcar Meneses Viveros on 13/03/19.
//
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv) {
    int rank, size;
    const int MAX_NUMBERS = 100;
    int numbers[MAX_NUMBERS];
    int number_amount;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
        fprintf(stderr, "Debe usar dos procesos\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        srand(time(NULL));
        number_amount = (rand() / (float)RAND_MAX) * MAX_NUMBERS;
        MPI_Send(numbers, number_amount, MPI_INT, 1, 10, MPI_COMM_WORLD);
        printf("0 sent %d numbers to 1\n", number_amount);
    } else if (rank == 1) {
        MPI_Recv(numbers, MAX_NUMBERS, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &number_amount);
        printf("1 received %d numbers from 0. Message source = %d, tag = %d\n",
               number_amount, status.MPI_SOURCE, status.MPI_TAG);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

