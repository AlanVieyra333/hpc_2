#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#define M 15000
#define N 15000

double a[N], *b, c[M];

int main()
{
    clock_t t_ini, t_fin;
    int i, j;
    double t1, t2;
    double sum;

    b = (double *)malloc(M * N * sizeof(double));
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            b[i * N + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (i = 0; i < N; i++)
        a[i] = i + 1;

    t_ini = clock();
    #pragma omp parallel shared(a,b,c) private(sum)
    {

        #pragma omp for private(j)
        for (i = 0; i < M; i++)
        {
            sum = 0.0;

            for (j = 0; j < N; j++)
            {
                sum += b[i * N + j] * a[j];
            }

            #pragma omp critical
            c[i]=sum;
        }
    }

    t_fin = clock();

//    printf("C = %lf\n", c[i]);

    for (i = 0; i < M; i++)
    {
        //printf("C = %lf ", c[i]);
    }

    //printf("\n");

    double secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    printf("%lf milisegundos\n", secs * 1000.0);

    free(b);

    return 0;
}
