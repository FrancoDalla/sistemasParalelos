#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MAX_SIZE 2000
#define COORDINATOR 0

int main(int argc, char* argv[]){
	int i, numProcs, rank, size, strip_size, local_sum=0, sum=0;
	int array[MAX_SIZE];
	MPI_Status status;

	size = atoi(argv[1]);
	size = (size < MAX_SIZE ? size : MAX_SIZE);

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == COORDINATOR)
		for(i=0; i<size; i++)
			array[i] = i+1;
	strip_size = size / numProcs;

	if(rank == COORDINATOR){
		for(i=1; i< numProcs; i++)
			MPI_Send(array+i*strip_size, strip_size, MPI_INT, i, 0, MPI_COMM_WORLD);
	} else
		MPI_Recv(array, strip_size, MPI_INT, COORDINATOR, 0, MPI_COMM_WORLD, &status);

	for(i=0; i<strip_size; i++)
		local_sum += array[i];

	if(rank == COORDINATOR) {
		sum = local_sum;
		for(i = 1; i<numProcs; i++){
			MPI_Recv(&local_sum, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			sum += local_sum;
		}
	}else
		MPI_Send(&local_sum, 1, MPI_INT, COORDINATOR, 1, MPI_COMM_WORLD);

	MPI_Finalize();

	if(rank == COORDINATOR)
		printf("SUMA = %d\n", sum);

	return 0;
}
