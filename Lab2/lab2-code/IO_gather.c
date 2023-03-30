#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

int main( argc, argv )
int  argc;
char **argv;
{
  int rank, size, i, ierr;
  char message[30];
  char *allmessages;
  int root = 0;
  
  MPI_Init( &argc, &argv );                 /* MPI functions return */
  MPI_Comm_size( MPI_COMM_WORLD, &size );   /* an error value that  */
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );   /* you should check!!   */
  
  sprintf(message, "hello world from proc %.2d", rank);

  if (rank == root) allmessages = (char *) malloc(size*30);

  MPI_Gather(message, 30, MPI_CHAR, 
             allmessages, 30, MPI_CHAR, root, MPI_COMM_WORLD);

  if (rank == root)
    for(i=0;i<size;i++)
      printf("Proc %.2d: %s\n", i, allmessages+30*i);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
