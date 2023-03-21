/**********************************************************************
 * A simple "hello world" program for MPI/C
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);               /* Initialize MPI               */

  printf("Hello World!\n");             /* Print a message              */

  MPI_Finalize();                       /* Shut down and clean up MPI   */

  return 0;
}
