//
//  pi_mpi.c
//  
//
//  Created by Alan Vieyra on 22/02/20.
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define NE  100000000

class Rand {
private:
    unsigned long int seed;
public:
    Rand();
    Rand(unsigned long int seed);
    int next(void);
};

Rand::Rand() {
    this->seed = 1;
}

Rand::Rand(unsigned long int seed) {
    this->seed = seed;
}

int Rand::next(void) // RAND_MAX assumed to be 32767
{
    this->seed = this->seed * 1103515245 + 12345;
    return (unsigned int)(this->seed/65536) % RAND_MAX;
}

long getMicrotime(){
	struct timeval currentTime;
	gettimeofday(&currentTime, NULL);
	return currentTime.tv_sec * (int)1e6 + currentTime.tv_usec;
}

long calcula_buenos_pi(long n, int rank) {
    long i;
    long eb=0;
    double x, y;

    Rand custom_rand(getMicrotime() ^ rank);
    
    srand(time(NULL));
    for (i=0; i<n; i++) {
        x = (double)custom_rand.next()/(double)RAND_MAX;
        y = (double)custom_rand.next()/(double)RAND_MAX;
        if ((x*x + y*y) <= 1.0) eb++;
    }
    
    return eb;
}

int main(int argc, char **argv) {
    int size,  rank;
    long eb, teb;
    int rc;
    double mi_pi;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    printf("La tarea %d se inicia\n", rank);
    
    eb = calcula_buenos_pi(NE, rank);
    
    rc = MPI_Reduce(&eb, &teb, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank==0) {
        mi_pi = 4.0*((double)teb/((double)size*NE));
        printf("Aproximado de PI = %lf\n", mi_pi);
        printf("El de math.h es %lf\n", M_PI);
    }
    
    MPI_Finalize();
    return 0;
}
