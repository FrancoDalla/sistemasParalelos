#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0){
        int vector[10] = {0, 1 , 2, 3, 4, 5, 6, 7, 8, 9};
        MPI_Ssend(vector+5, 5, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    else{
        int *vector = (int *) malloc(sizeof(int) * 5);
        MPI_Recv(vector, sizeof(int) * 5, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int i = 0; i < 5; i++){
            printf("%d ",vector[i]);
        }
        printf("\n");
    }

    MPI_Finalize();
    return 0;
}
