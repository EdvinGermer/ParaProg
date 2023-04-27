#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void comm_splitR(int local_max, int local_min, int level, MPI_Comm communicator);

int main(int argc, char **argv) {

	// Parallel environment
	MPI_Init(&argc, & argv);
	int rank, num_proc;
	int local_max,local_min,level;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	
	// Input data on local PEs
	local_max = (rank+1)*100;
	local_min = rank+1;
	level = 0;
	int try = 1;
	while(1){
	      if(num_proc == try){
	      break;
	      }
	      try = try*2;
	      level = level + 1;
	 }
        printf("P:%d (before splitting) Level: %d, max= %d, min= %d \n", rank,level,local_max,local_min);
	comm_splitR(local_max,local_min,level, MPI_COMM_WORLD);
	
        MPI_Finalize();
	return 0;
	}

void comm_splitR(int local_max, int local_min, int level, MPI_Comm communicator) {
	int rank, num_proc;
	int global_max = 0, global_min = 0;
	
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &num_proc);
	// Base case: 1 processors in communicator
	if (1 == num_proc) {
		return;

	} else {
                MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, communicator);	        
		MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, communicator);
                printf("P:%d Level: %d, max= %d, min= %d \n", rank,level,global_max,global_min);		
		// Split communicator
		MPI_Comm halved_communicator;
		int color = (rank >= num_proc/2);
		MPI_Comm_split(communicator, color, rank, &halved_communicator);
                level = level - 1;
		
		// Call recursively
		comm_splitR(local_max, local_min, level, halved_communicator);
		MPI_Comm_free(&halved_communicator);
		return;
	}
}
