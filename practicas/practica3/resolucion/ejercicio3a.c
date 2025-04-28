#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include <sys/time.h>

double dwalltime();

int main(int argc, char *argv[]){
    double *A, *B, *C;
    int i, j, k, N;
    int check = 1;
    double timetick;

    N = atoi(argv[1]);
    int numThreads = atoi(argv[2]);
    omp_set_num_threads(numThreads);

    A = (double*) malloc(sizeof(double)*N*N);
    B = (double*) malloc(sizeof(double)*N*N);
    C = (double*) malloc(sizeof(double)*N*N);

    /*Inicializacion */
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            A[i*N+j] = 1;
            B[i+j*N] = 1;
        }
    }

    timetick = dwalltime();

    /*multiplicacion */
    #pragma omp parallel for
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            C[i*N+j] = 0;
            for(k = 0; k < N; k++){
                C[i*N+j] = C[i*N+j] + A[i*N+k] * B[k+j*N];
            }
        }
    }

    printf("Tiempo en segundos %f\n",
        dwalltime() - timetick);


}



double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}
