/************************************************************************
 * This file has been written as a sample solution to an exercise in a 
 * course given at the Edinburgh Parallel Computing Centre. It is made
 * freely available with the understanding that every copy of this file
 * must include this header and that EPCC takes no responsibility for
 * the use of the enclosed teaching material.
 *
 * Author:      Joel Malard     
 *
 * Contact:     epcc-tec@epcc.ed.ac.uk
 *
 * Purpose:     A program to try out non-blocking point-to-point 
 *              communications.
 *
 * Contents:    C source code.
 *
 ************************************************************************/

#include        <stdio.h>
#include        <mpi.h>
#define to_right 201
#define to_left  102

void main (int argc, char *argv[])
{
  int ierror, rank, my_rank, size;
  int right, left;
  int other, sum, i;
  
  MPI_Status  send_status;
  MPI_Status  recv_status;
  MPI_Request request;
  
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  right = my_rank + 1;
  if (right == size) right = 0;
  
  left = my_rank - 1;
  if (left == -1) left = size-1;
  
  sum = 0;
  rank = my_rank;
  
  for( i = 0; i < size; i++) {
    
    MPI_Issend(&rank, 1, MPI_INT, right, to_right, MPI_COMM_WORLD, &request);
    
    MPI_Recv(&other, 1, MPI_INT, left, to_right, MPI_COMM_WORLD, &recv_status);
    
    MPI_Wait(&request, &send_status);
    
    sum = sum + other;
    rank = other;
    
  }
  
  printf ("PE%d:\tSum = %d\n", rank, sum);
  
  MPI_Finalize();

}
