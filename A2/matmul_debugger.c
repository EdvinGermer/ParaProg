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
    double *A_local, *B_local, *C_local;
    
    /* CHECK ARGUMENT COUNT */
    if (argc != 3)
    {
        if (rank == 0)
            printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return -1;
    }

    /* INITIALIZE MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

        // Check matrices were read correctly
        printf("Matrix A:\n");
        print_matrix(A, n,n);
        printf("\n\nMatrix B:\n");
        print_matrix(B, n,n);
        printf("__________________________\n");
    }



    /* DISTRIBUTE DATA */
    double start = MPI_Wtime(); // Record the start time

    // Broadcast the value of n to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //printf("Rank %d out of %d\n",rank,size);

    // Determine block size
    int m = n/size;

    // Create local blocks
    A_local = (double *)malloc(m * n * sizeof(double));
    B_local = (double *)malloc(n * m * sizeof(double));
    C_local = (double *)malloc(m * n * sizeof(double));

    // Send block of rows of A from rank 0 to A_local on each rank
    MPI_Scatter(A, n * m, MPI_DOUBLE, A_local, n * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //printf("Rank %d, matrix A:\n",rank);
    //print_matrix(A_local, m,n);

    // Send block of cols of B from rank 0 to B_local on each rank
    if (rank == 0)
    {
        // Copy data from big B matrix to B_local
        for (int j = 0; j < n; j++)
            memcpy(&B_local[j * m], &B[j * n], m * sizeof(double));

        // Send data to other processes
        for (int i = 1; i < size; i++) // Other processes
            for (int j = 0; j < n; j++) // Element index
                MPI_Send(&B[j * n + i * m], m, MPI_DOUBLE, i, j, MPI_COMM_WORLD);
    }
    else // Other processes recieves data from 0
        for (int j = 0; j < n; j++)
            MPI_Recv(&B_local[j * m], m, MPI_DOUBLE, 0, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //printf("Rank %d, matrix B:\n",rank);
    //print_matrix(B_local, n,m);

    if (rank==0)
        printf("__________________________\n");
    


    /* PERFORM MATRIX MULTIPLICATION */
    // A_local has size (m,n) = m*n elements
    // B_local has size (n,m) = m*n elements
    // C_local has size (m,n) = size*m*n elements
    // Each stage calculates (m,m) elements in local

    for (int stage=0; stage<size; stage++) // Iterate over all stages (size-1 stages)
    {   
        // Matrix multiplication for current stage
        double *block = (double *)malloc(m * m * sizeof(double));
        for (int i = 0; i < m; i++) // Rows
            for (int j = 0; j < m; j++) // Cols
                {
                    double sum = 0.0;
                    for (int k = 0; k < n; k++) // Sum up
                    {
                        sum += A_local[i * n + k] * B_local[k * m + j];
                    }
                    block[i*m+j] = sum;
                }

        // Print block matrix
        printf("Rank %d, block matrix:\n",rank);
        print_matrix(block,m,m);

        // Copy block into C_local
        int start_col = (m * (rank+stage)) % n;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                C_local[i * n + (start_col + j)] = block[i * m + j];
            }
        }
        free(block);

        // Print C_local matrix
        printf("Rank %d, c_local matrix:\n",rank);
        print_matrix(C_local,m,n);

        // Shift B_local one step
        double *B_temp = (double *)malloc(n * m * sizeof(double));
        int prev = (rank - 1 + size) % size; // recieve from previous
        int next = (rank + 1) % size; // Send to next
        MPI_Sendrecv(B_local, m * n, MPI_DOUBLE, next, 1, B_temp, m * n, MPI_DOUBLE, prev, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


        // Overwrite old B_local with temp and free memory
        memcpy(B_local, B_temp, m * n * sizeof(double));
        free(B_temp);
    }
    
    /* GATHER RESULTS BACK TO RANK 0 */
    if (rank==0)
        C = (double *)malloc(n * n * sizeof(double));
    MPI_Gather(C_local, m * n, MPI_DOUBLE, C, m * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* PRINT FINAL CALCULATION */
    if (rank == 0)
    {
        printf("__________________________\n");
        printf("Result Matrix C:\n");
        print_matrix(C, n, n);
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

    MPI_Finalize();
    return 0;
}