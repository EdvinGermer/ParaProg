#include<mpi.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

// mpirun -np 4 ./cart_shift_NSEW

int main(int argc, char *argv[])
{
    int errs = 0, i;
    int me, numprocs, source, north, east, south, west;
    int nw, ne, sw, se;
    int ndims, dims[2], periods[2], reorder; 
    int coord[2], coord_ne[2], coord_sw[2], coord_nw[2], coord_se[2];
    int mx = 10;
    int m[10];
    int n[10];
    
    MPI_Comm comm2D;
    MPI_Status status;
    
/* Beginning of the program */
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    
    ndims   = 2;
    dims[0] = sqrt(numprocs); dims[1] = sqrt(numprocs);
    periods[0] = 1; periods[1] = 1; reorder=0; 

    MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,periods,reorder,&comm2D);
        
    MPI_Cart_shift(comm2D,0,1,&source,&east);
    MPI_Cart_shift(comm2D,0,-1,&source,&west);
    MPI_Cart_shift(comm2D,1,1,&source,&north);
    MPI_Cart_shift(comm2D,1,-1,&source,&south);
    printf("For %d north-south-me-east-west:  %d <-- %d <-- %d -->%d -->%d \n", me,north,south,me,east,west); 

/*    nw=west;
    printf("P:%d value nw is set to: %d \n", me,nw);
    printf("P:%d source: %d \n", me,source);
    MPI_Send(&west,1,MPI_INT,south,222,comm2D);
    MPI_Recv(&nw,1,MPI_INT,north,222,comm2D,&status);
    printf("P:%d My NW neighbor is: %d \n", me,nw); 
*/      
    printf("P:%d Source array: ", me);
    for(i=0; i<mx; i++) {
        m[i] = me;
	n[i] = me+100;
        printf("  %d", m[i]);
    }
      printf("\n");

/*   MPI_Sendrecv_replace(m, mx, MPI_INT, north, 123, south, 124, comm2D, &status); 
   
    printf("P:%d Current  array:", me);
    for(i=0; i<mx; i++) {
            printf("  %d", m[i]);
    }
      printf("\n");*/
      

    // COORDS
    MPI_Cart_coords(comm2D, me, 2, coord);
    printf("P:%d coord: %d  %d \n", me,coord[0],coord[1]); 
    
    // FIND SW AND NE
    coord_ne[0]=coord[0]+1; coord_ne[1]=coord[1]+1;
    MPI_Cart_rank(comm2D, coord_ne, &ne); 
    coord_sw[0]=coord[0]-1; coord_sw[1]=coord[1]-1;
    MPI_Cart_rank(comm2D, coord_sw, &sw); 
    printf("P:%d  sw-me-ne:  %d <-- %d -->%d \n", me,sw,me,ne); 

    // FIND SE AND NW
    coord_nw[0]=coord[0]-1; coord_nw[1]=coord[1]+1;
    MPI_Cart_rank(comm2D, coord_nw, &nw); 
    coord_se[0]=coord[0]+1; coord_se[1]=coord[1]-1;
    MPI_Cart_rank(comm2D, coord_se, &se); 
    printf("P:%d  se-me-nw:  %d <-- %d -->%d \n", me,se,me,nw); 

/*    MPI_Sendrecv_replace(m, mx, MPI_INT, ne, 123, sw, 124, comm2D, &status);*/
    MPI_Sendrecv(m, mx, MPI_INT, ne, 123, n, mx, MPI_INT, sw, 123, comm2D, &status);
    printf("P:%d Received  array:", me);
    for(i=0; i<mx; i++) {
            printf("  %d", n[i]);
    }
     printf("\n");
    
      
MPI_Comm_free( &comm2D );


MPI_Finalize();

return 0;

}
