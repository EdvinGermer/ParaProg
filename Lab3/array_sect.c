/*
 * On src:   m(0:9,0:9) 
 *           send section m(2:6,4:9)
 *
 * On dest   m(0:11,0:7)
 *           recv section m(7:11,2:7)
 *
 * NOTE: this code emulates Fortran style indexing
 */
#include  <stdio.h>
#include  "mpi.h"

#define IND(i,j)  j*mx+i

void main( int argc, char *argv[] )
{
  int *m;
  int mx, my, nx, ny;
  int ierr, rank, size;
  int src = 1, dest = 0;
  MPI_Datatype section;
  MPI_Status status;
  int i,j;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank==src) {

    mx = 10; my = 10;
    m = (int *) malloc(mx*my*sizeof(int));

    printf("Source array:\n");
    for(j=0; j<my; j++) {
      for(i=0; i<mx; i++){
	m[IND(i,j)] = IND(i,j);
	printf("  %2.2d", m[IND(i,j)]);
      }
      printf("\n");
    }
    printf("\n");

    nx = 6-2+1; ny = 9-4+1;  /* m(2:6,4:9) */
  }
  else{

     mx = 12; my = 8;
     m = (int *) malloc(mx*my*sizeof(int));
     for(j=0; j<mx*my; j++) m[j] = 0;

     nx = 11-7+1; ny = 7-2+1;  /* m(7:11,2:7) */
  }

  MPI_Type_vector(ny, nx, mx, MPI_INT, &section);

  MPI_Type_commit(&section);
  
  if (rank==src)
    MPI_Send( m+IND(2,4), 1, section, dest, 7, MPI_COMM_WORLD);
  else
    {
      MPI_Recv( m+IND(7,2) , 1, section, src, 7, MPI_COMM_WORLD, &status);

      printf("Target array:\n");
      for(j=0; j<my; j++) {
	for(i=0; i<mx; i++) printf("  %2.2d", m[IND(i,j)]);
	printf("\n");
      }
      printf("\n");
      
    }

  MPI_Finalize();

}


