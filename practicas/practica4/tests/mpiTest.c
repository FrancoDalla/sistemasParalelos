#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
	int cantidad, identificador;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &cantidad);
	MPI_Comm_rank(MPI_COMM_WORLD, &identificador);
	
	if(identificador == cantidad - 1)
	{
		printf("HOLA SOY EL ULTIMO RANGOOOOO\n", identificador);
	}
	else
	{
		printf("HOLA MUNDO SOY %d de %d \n",
		identificador, cantidad);
	}
	printf("%d\n",cantidad);
	MPI_Finalize();
	printf("%d\n",cantidad);
	return 0;
}
