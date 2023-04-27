 /**********************************************************************
 * Derived datatypes in MPI/C, transpose a 2D array
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
	// Check arguments
	if (3 != argc) {
		printf("Usage: ./datatypes_transpose input1 input2, where input1 is the number of rows and input2 is the number of columns in the local array.\n");
		return -3;
	}
	char *nnx = argv[1];
	char *nny = argv[2];
	
  int rank,size,row,col,count,blocklen,stride, nx,ny;
  int errs = 0;
  int source, left, right;
  int ndims, dims[1], periods[1], reorder;   
  double *A, *B;
  nx = atoi(nnx); ny = atoi(nny);
 
  MPI_Status status;
  MPI_Datatype rowtype;
  MPI_Datatype coltype;
  MPI_Comm comm1D;
  MPI_Request reqs[nx];
  MPI_Status stats[nx];

  MPI_Init(&argc, & argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
  
  ndims   = 1;
  dims[0] = size; 
  periods[0] = 1; reorder=0; 
  MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&comm1D);  
  MPI_Cart_shift(comm1D,0,1,&source,&right);
  MPI_Cart_shift(comm1D,0,-1,&source,&left);
  
  A=(double *)malloc(nx*ny*sizeof(double));
  B=(double *)malloc(nx*ny*sizeof(double));
   
  for (row=0; row<nx;row++)
  {
    for (col=0; col<ny;col++)
    {
      A[row*ny+col]=(double)row*10+(double)col+1+rank;
      B[row*ny+col]=0.0;
/*      printf("%d: In: A=%f,  B=%f\n",rank,A[row*ny+col],B[row*ny+col]); */
    }
  }  

  count=nx; blocklen=1; stride=ny;
  MPI_Type_vector(count,blocklen,stride,MPI_DOUBLE,&coltype);
  MPI_Type_commit(&coltype);
  MPI_Type_contiguous(ny,MPI_DOUBLE,&rowtype);
  MPI_Type_commit(&rowtype);
  
 
   for (row=0; row<nx; row++)
   {
 /*      MPI_Irecv(&B[row*ny+1], 1,coltype, left,  right+row, comm1D, &reqs[row]);*/
        MPI_Send(&A[row], 1, coltype, right, right+row, comm1D);
        MPI_Recv(&B[row*ny], 1, rowtype, left,  rank+row, comm1D, &stats[row]);
  }
/*   MPI_Waitall(5, reqs, stats); */
 
 
   
   for (row=0; row<nx; row++)
   {
       for (col=0; col<ny; col++)
       {
      printf("%d: Out: %d: A=%f,  B=%f\n",rank,row,A[row*ny+col],B[row*ny+col]);
      }
   }
   
  free(A); free(B);
  MPI_Type_free(&rowtype);
  MPI_Type_free(&coltype);
  MPI_Comm_free(&comm1D);

  MPI_Finalize(); 

  return 0;
}
