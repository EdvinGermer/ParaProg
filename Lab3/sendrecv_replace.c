#include <mpi.h>
#include <math.h>
#include <stdio.h>
 
int main(int argc, char *argv[])
{
    int errs = 0;
    int me, numprocs, source, north, east, south, west;
    int ndims, dims[2], periods[2], reorder; 
    int mybuf[2];
    MPI_Comm comm2D;
    MPI_Status status;
    
/* Beginning of the program */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
        mybuf[0]=me; mybuf[1]=me+1; 
        ndims   = 2;
        dims[0] = sqrt(numprocs); dims[1] = sqrt(numprocs);
        periods[0] = 1; periods[1] = 1; 
        reorder=0; 

      MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&comm2D);
        
      MPI_Cart_shift(comm2D,0,1,&source,&east);
      MPI_Cart_shift(comm2D,0,-1,&source,&west);
      MPI_Cart_shift(comm2D,1,1,&source,&north);
      MPI_Cart_shift(comm2D,1,-1,&source,&south);
      printf("For %d north-west-me-east-south:  %d <-- %d <-- %d -->%d  -->%d \n", me,north,west,me,east,south);

      printf("Before: for %d: mybuf[%d, %d] \n", me, mybuf[0], mybuf[1]);
      MPI_Sendrecv_replace(mybuf, 2, MPI_INT, west, 99, east, 99, comm2D, &status);
      printf("After: for %d: mybuf[%d, %d] \n", me, mybuf[0], mybuf[1]);

      MPI_Comm_free( &comm2D );
      
      MPI_Finalize();
      return 0;
}
