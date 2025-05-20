/*
** Sending simple, point-to-point messages.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"

#define MASTER 0

int main(int argc, char* argv[])
{
  int myrank;
  int size;
  int dest;              /* destination rank for message */
  int source;            /* source rank of a message */
  int tag = 0;           /* scope for adding extra information to a message */
  MPI_Status status;     /* struct used by MPI_Recv */
  char message[BUFSIZ];

  /* MPI_Init returns once it has started up processes */
  MPI_Init( &argc, &argv );

  /* size and rank will become ubiquitous */
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

  if(myrank == MASTER){
      source = size - 1;
      dest = myrank + 1;
      char send_message[BUFSIZ];
      sprintf(send_message, "%d", myrank);

      MPI_Send(send_message, strlen(message) + 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
      MPI_Recv(message, BUFSIZ, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
  }
  else{
      if(myrank == size - 1){
          dest = MASTER;

      }
      else{
          dest = myrank + 1;
      }

      source = myrank - 1;

      char send_message[BUFSIZ];
      sprintf(send_message, "%d", myrank);

      MPI_Recv(message, BUFSIZ, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
      MPI_Send(send_message, strlen(message) + 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
  }

  printf("Soy el proceso %d y se comunico conmigo el proceso %s\n", myrank, message);

  MPI_Finalize();

  /* and exit the program */
  return EXIT_SUCCESS;
}
