#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MASTER 0

MPI_Status status;
MPI_Request req;

int *a;
int n;

void jefe_creacion_datos();
void jefe_envio_datos();
void delegado_recepcion();
void delegado_trabajo();

int size;
int myrank;
int strip_size;

int main(int argc, char *argv[]){
    int dest;
    int source;
    int tag = 0;
    n = atoi(argv[1]);


    char message[BUFSIZ];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(myrank == 0){
        jefe_creacion_datos();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(myrank == 0){
        jefe_envio_datos();
    }
    else{
        delegado_recepcion();
    }

    delegado_trabajo();

    MPI_Finalize();
    return 0;
}

void jefe_creacion_datos(){
    strip_size = n / size;
    a = (int *) malloc(sizeof(int) * n * n);
    for(int i = 0; i <n ; i++){
        for(int j= 0; j <n; j++){
            a[i*n+j] = 1;
        }
    }
}

void jefe_envio_datos(){
    for(int i = 1; i<size; i++){
        MPI_Isend(a + i * strip_size * n, strip_size * n, MPI_INT, i, 0, MPI_COMM_WORLD, &req);
    }
}

void delegado_recepcion(){
    MPI_Irecv(a, strip_size * n, MPI_INT, 0, 0,MPI_COMM_WORLD, &req);
    MPI_Wait(&req, &status);
}

void delegado_trabajo(){
    for(int i = 0; i< strip_size; i++){
        printf("%d, ", a[i]);
        printf("\n");
    }
}
