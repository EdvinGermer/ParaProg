/**********************************************************************
 * A simple "hello world" program for MPI/C
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {

  int size,rank;
  
  MPI_Init(&argc, &argv);               /* Initialize MPI               */

  printf("Hello World!\n");             /* Print a message              */

  // ADDED
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Write current rank to "rank"
  MPI_Comm_size(MPI_COMM_WORLD, &size);  // Write total # processes to "size"
  
  printf("  Rank: %d\n",rank); 
  if (rank==0)
    printf("  Total number of processes: %d\n", size); 

  MPI_Finalize();                       /* Shut down and clean up MPI   */

  return 0;
}
