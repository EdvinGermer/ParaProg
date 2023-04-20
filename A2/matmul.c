#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include<unistd.h>

void print_matrix(double *matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.2lf ", matrix[i * cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// mpirun --bind-to none -n 2 ./matmul ./input4.txt output.txt

int main(int argc, char *argv[])
{
    /* INITIALIZE PARAMETERS */
    int rank,size,n;
    double *A, *B, *C;
    double *A_local, *B_local, *B_temp, *C_local;

    /* INITIALIZE MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* CHECK ARGUMENT COUNT */
    if (argc != 3)
    {
        if (rank == 0)
            printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return -1;
    }

    /* READ DATA TO RANK 0 */
    if (rank == 0)
    {
        FILE *input = fopen(argv[1], "r");

        // Read the size of the matrices from the input file
        if (fscanf(input, "%d", &n) != 1)
        {
            printf("Error when reading n from input file");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Allocate memory for matrix A and read its values from the input file
        A = (double *)malloc(n * n * sizeof(double));
        for (int i = 0; i < n * n; i++)
        {
            if (fscanf(input, "%lf", &A[i]) != 1)
            {
                printf("Error when reading A from input file");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // Allocate memory for matrix B and read its values from the input file
        B = (double *)malloc(n * n * sizeof(double));
        for (int i = 0; i < n * n; i++)
        {
            if (fscanf(input, "%lf", &B[i]) != 1)
            {
                printf("Error when reading B from input file");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
        fclose(input); // Close the input file
    }

    /* DISTRIBUTE DATA */
    double start = MPI_Wtime(); // Record the start time

    // Broadcast the value of n to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Determine block size
    int m = n/size;

    // Create local blocks
    A_local = (double *)malloc(m * n * sizeof(double));
    B_local = (double *)malloc(n * m * sizeof(double));
    B_temp = (double *)malloc(n * m * sizeof(double));
    C_local = (double *)malloc(m * n * sizeof(double));

    // Send block of rows of A from rank 0 to A_local on each rank (figure 2 in report)
    MPI_Scatter(A, m*n, MPI_DOUBLE, A_local, m*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Send block of cols of B from rank 0 to B_local on each rank  (figure 2 in report)
    if (rank == 0)
    {
        // Copy data from big B matrix to B_local
        for (int j = 0; j < n; j++)
            memcpy(&B_local[j * m], &B[j * n], m * sizeof(double)); //m elements from b[j*n]

        // Send data to other processes
        for (int i = 1; i < size; i++) // Other processes
            for (int j = 0; j < n; j++) // Iterate over rows
                {
                    // j*n = row index
                    // i*m = where column block start on that row
                    // send m elements per MPI_send (column thickness)
                    MPI_Send(&B[j * n + i * m], m, MPI_DOUBLE, i, j, MPI_COMM_WORLD); // from 0 to i
                } 
    }
    else // Other processes recieves data from 0
        for (int j = 0; j < n; j++)
            MPI_Recv(&B_local[j * m], m, MPI_DOUBLE, 0, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    /* PERFORM MATRIX MULTIPLICATION */
    for (int stage=0; stage<size; stage++) // Iterate over all stages
    {   
        // Index where to write each block in the local C_matrix
        int start_col = (m * (rank+stage)) % n;  // Start with diagonal then step right in C
        
        // Matrix multiplication for current stage
        // Block of rows from A matrix multiplied with block columns from B Matrix
        for (int i = 0; i < m; i++) // Rows
            for (int j = 0; j < m; j++) // Cols
                {
                    double sum = 0.0;
                    for (int k = 0; k < n; k++) // Sum up
                    {
                        sum += A_local[i * n + k] * B_local[k * m + j];
                    }
                    C_local[i * n + (start_col + j)] = sum;  // i*n = select row,   start_col+j = block width is start_col to start_col+m
                }

        // Shift B_local one step (to the left)
        int dest = (rank - 1 + size) % size; // send to previous
        int source = (rank + 1) % size; // recieve from next
        MPI_Sendrecv(B_local, m * n, MPI_DOUBLE, dest, 1, B_temp, m * n, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Overwrite old B_local with temp and free memory 
        memcpy(B_local, B_temp, m * n * sizeof(double));
    }
    
    /* GATHER RESULTS BACK TO RANK 0 */
    if (rank==0)
        C = (double *)malloc(n * n * sizeof(double));
    MPI_Gather(C_local, m * n, MPI_DOUBLE, C, m * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* LOCAL TIME */
	double local_execution_time = MPI_Wtime() - start;

    /* WRITE RESULTS TO OUTPUT */
    if (rank==0)
    {
        FILE *output = fopen(argv[2], "w");
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                fprintf(output, "%.6lf ", C[i * n + j]);
            }
            fprintf(output, "\n");
        }
        fclose(output);
    }

    /* FIND LONGEST TIME */
	double timings[size];
	MPI_Gather(&local_execution_time, 1, MPI_DOUBLE, timings, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		double max = 0;
		for (int i=0; i<size;i++)
		{
			double time = timings[i];
			if (time>max)
				max = time;
		}
		printf("%f\n", max);
	}
    
    /* FREE MEMORY */
    if (rank==0)
    {
        free(A);
        free(B);
        free(C);
    }
    
    free(A_local);
    free(B_local);
    free(C_local);
    free(B_temp);

    MPI_Finalize();
    return 0;
}