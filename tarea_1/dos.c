//
//  status2.c
//  
//
//  Created by Amilcar Meneses Viveros on 13/03/19.
//
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv) {
    int size, rank;
    int number_amount;
    const int MAX_NUMBERS = 100;
    int numbers[MAX_NUMBERS];
    MPI_Status status;
    int *number_buf;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
        fprintf(stderr, "Must use two processes for this example\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        srand(time(NULL));
        number_amount = (rand() / (float)RAND_MAX) * MAX_NUMBERS;
        MPI_Send(numbers, number_amount, MPI_INT, 1, 10, MPI_COMM_WORLD);
        printf("0 sent %d numbers to 1\n", number_amount);
    } else if (rank == 1) {
        MPI_Probe(0, 10, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &number_amount);
        number_buf = (int *)malloc(sizeof(int) * number_amount);
        MPI_Recv(number_buf, number_amount, MPI_INT, 0, 10, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        printf("1 dynamically received %d numbers from 0.\n",
               number_amount);
        free(number_buf);
    }
    MPI_Finalize();
}

