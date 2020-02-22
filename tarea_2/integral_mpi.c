//
//  inrtegral_mpi.c
//  
//
//  Created by Amilcar Meneses Viveros on 19/02/20.
//
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NE   1000000

#define XMIN   -M_PI
#define XMAX    M_PI 

double sin(double);
double cos(double);
double exp(double);

double func(double x) {
    return sin(x)*cos((x*x)*sin(x*x));
}

double integral(int xmin, int xmax, int n, double (*f)(double)) {
    int i;
    double tmp=0;
    double x, dx = (xmax-xmin)/(double)n;
    
    x = xmin;
    for (i=0; i<n; i++) {
        tmp += (*f)(x)*dx;
        x+=dx;
    }
    return tmp;
}

double gaussiana(double x) {
    return exp(-x*x);
}

int main(int argc, char **argv) {
    int size,  rank;
    int rc;
    double xmin, xmax, dx;
    double inttmp, tint;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    printf("La tarea %d se inicia\n", rank);
    
    dx = (XMAX-XMIN)/size;
    
    xmin = XMIN + rank*dx;
    xmax = xmin+dx;
    inttmp = integral(xmin, xmax, NE/size, func);
    
    rc = MPI_Reduce(&inttmp, &tint, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank==0) {
        printf("La integral es = %lf\n", tint);
    }
    
    MPI_Finalize();
    return 0;
}
