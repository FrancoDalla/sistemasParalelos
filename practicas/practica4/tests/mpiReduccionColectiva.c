#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MAX_SIZE 2000
#define COORDINATOR 0

int main(int argc, char* argv[]){
    int i, num_procs, rank, size, strip_size, local_sum=0, sum=0;
    int array[MAX_SIZE];
    MPI_Status status;

    size = atoi(argv[1]);
    size = (size < MAX_SIZE ? size : MAX_SIZE);

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == COORDINATOR)
        for(i=0; i<size; i++)
            array[i] = i+1;

    strip_size = size / num_procs;

    MPI_Scatter(array, strip_size, MPI_INT, array, strip_size, MPI_INT, COORDINATOR, MPI_COMM_WORLD);

    for(i=0; i<strip_size; i++)
        local_sum += array[i];

    MPI_Reduce(&local_sum, &sum, 1, MPI_INT, MPI_SUM, COORDINATOR, MPI_COMM_WORLD);

    if(rank==COORDINATOR)
        printf("SUM = %d\n",sum);

    MPI_Finalize();

    return 0;
}
