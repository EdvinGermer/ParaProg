#include <mpi.h>
#include <stdio.h>

#define ROOT 0

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int rank, number;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (ROOT == rank) {
		printf("Enter a number! Each process will print the sum of the number and its own rank.\n");
		scanf("%d", &number);
		MPI_Bcast(&number, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	printf("P%d: %d\n", rank, rank+number);
	MPI_Finalize();
	return 0;
}
