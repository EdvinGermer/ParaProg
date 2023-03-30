#include <mpi.h>
#include <stdlib.h>

#define N 3

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// One-dimensional array
	double array1D[N];
	for (int i=0; i<N; i++) {
		array1D[i] = rank + (double)rand()/RAND_MAX;
	}

	// Statically allocated 3-dimensional array
	double array3D_static[N][N][N];
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			for (int k=0; k<N; k++) {
				 array3D_static[i][j][k] = rank + (double)rand()/RAND_MAX;
			 }
		}
	}

	// Dynamically allocated 3-dimensional array
	double ***array3D_dynamic = malloc(N*sizeof(double **));
	for (int i=0; i<N; i++) {
		array3D_dynamic[i] = malloc(N*sizeof(double *));
		for (int j=0; j<N; j++) {
			array3D_dynamic[i][j] = malloc(N*sizeof(double));
			for (int k=0; k<N; k++) {
				array3D_dynamic[i][j][k] = rank + (double)rand()/RAND_MAX;
			}
		}
	}

	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			free(array3D_dynamic[i][j]);
		}
		free(array3D_dynamic[i]);
	}
	free(array3D_dynamic);

	MPI_Finalize();
	return 0;
}
