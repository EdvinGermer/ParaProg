/**********************************************************************
 * Point-to-point communication using MPI
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int rank, size;
  double a, b;
  MPI_Status status;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
  
  a = 100.0 + (double) rank;  /* Different a on different processors */

  /* Exchange variable a, notice the send-recv order */
  if (rank == 0) {
    MPI_Send(&a, 1, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD);
    MPI_Recv(&b, 1, MPI_DOUBLE, size-1, size-1, MPI_COMM_WORLD, &status);
    printf("Processor %d got %f from processor %d\n",rank, b,size-1);
  }
  else {
    for (int r=1;r<size-1;r++)
    {
      if (r==rank)
      {
        MPI_Recv(&b, 1, MPI_DOUBLE, r-1, r-1, MPI_COMM_WORLD, &status);
        MPI_Send(&a, 1, MPI_DOUBLE, r+1, r, MPI_COMM_WORLD);
        printf("Processor %d got %f from processor %d\n",r, b,r-1);
      }
    }

    // Last send/recieve
    if (rank==size-1)
    {
      MPI_Recv(&b, 1, MPI_DOUBLE, size-2, rank-1, MPI_COMM_WORLD, &status);
      MPI_Send(&a, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
      printf("Processor %d got %f from processor %d\n",rank, b,size-2);
    }
  }
  
  MPI_Finalize(); 

  return 0;
}
