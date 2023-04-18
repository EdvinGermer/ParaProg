#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include<unistd.h>

// Files @ /crex/proj/uppmax2023-2-13/nobackup/matmul_indata/input4.txt

// module load gcc openmpi
// make matmul
// ./matmul.sh

void print_block(double *block, int block_size, int rank) {
    printf("Rank %d:\n", rank);
    for (int i = 0; i < block_size; i++) {
        for (int j = 0; j < block_size; j++) {
            printf("%.2lf ", block[i * block_size + j]);
        }
        printf("\n");
    }
    printf("\n");
}


int main(int argc, char *argv[])
{
    /* INITIALIZE PARAMETERS */
    int rank,size;
    
    /* CHECK ARGUMENT COUNT */
    if (argc != 3)
    {
        if (rank == 0)
        {
            printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        }
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
    }







    MPI_Finalize();
    return 0;
}