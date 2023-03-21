/**********************************************************************
 * Point-to-point communication using MPI
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int rank, size,i;
  double a;
  MPI_Status status;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
  
  /* Processor 0 send to all others */
  if (rank == 0) {
    a=999.999;
    for (i=0;i<size;i++)
      MPI_Send(&a, 1, MPI_DOUBLE, i, 111, MPI_COMM_WORLD);
  } else {
    MPI_Recv(&a, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD, &status);
    printf("Processor %d got %f\n", rank,a);
  }

  MPI_Finalize(); 

  return 0;
}
