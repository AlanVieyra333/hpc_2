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

int main(int argc, char **argv) {
    int size, rank;
    int numbers_len = 30;
    int numbers[numbers_len];
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
        // Llenar array
        int i = 0;
        for(i =0; i<numbers_len; i++) {
            numbers[i] = i+1;
        }
        MPI_Send(numbers, numbers_len, MPI_INT, 1, 10, MPI_COMM_WORLD);
        printf("0 sent %d numbers to 1\n", numbers_len);
    } else if (rank == 1) {
        MPI_Probe(0, 10, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &numbers_len);
        number_buf = (int *)malloc(sizeof(int) * numbers_len);
        MPI_Recv(number_buf, numbers_len, MPI_INT, 0, 10, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        printf("1 dynamically received %d numbers from 0.\n",
               numbers_len);

        int i=0;
        for(i=0; i<numbers_len; i++){
            printf("%d ", number_buf[i]);
        }
        printf("\n");

        free(number_buf);
    }
    MPI_Finalize();
}
