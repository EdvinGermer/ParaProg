#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// module load gcc openmpi
// mpirun --bind-to none -n 2 ./exercise1

int main(int argc, char** argv)
{
    /* INITIALIZE PARAMETERS */
    int p, rank, k;
    int n = 6;      // big array size, where n is from the nxn matrix
    int m = 2;      // How many rows to send

    int *big_array, *local_array;

    /* INITIALIZE MPI*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    k = n / p;  // Local block width

    /* CREATE BIG ARRAY REPRESENTING nxn MATRIX */
    if (rank==0)
    {
        big_array = (int*)malloc(n*n * sizeof(int));
        for (int i=0; i<n*n; ++i)
            big_array[i] = i+1;
    }
    

    /* SCATTER BIG_ARRAY TO LOCAL_ARRAY */
    local_array = (int*)malloc(k*n * sizeof(int));
    MPI_Scatter(big_array, k*n, MPI_INT, local_array, k*n, MPI_INT, 0, MPI_COMM_WORLD);

    printf("RANK %d\n    ",rank);
    for (int i=0;i<k*n;i++)
        printf("%d, ",local_array[i]);
    printf("\n\n");
    

    
    

    /* DEFINE CONTIGUOUS DATATYPE FOR SENDING */
    MPI_Datatype contiguous;
    MPI_Type_contiguous(m * n, MPI_INT, &contiguous);  // m rows, n columns,   m<k
    MPI_Type_commit(&contiguous);


    /* SEND DATA */
    int right = (rank + 1) % p;
    int left = (rank - 1 + p) % p;

    int *temp = (int*)malloc(m*n * sizeof(int)); // temp array to recieve data from left neighbour
    int start_idx = k*n - m*n; // k*n = last element
    MPI_Sendrecv(&local_array[start_idx], 1, contiguous, right, 0, temp, 1, contiguous, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("RANK %d recieved:\n    ",rank);
    for (int i=0;i<m*n;i++)
        printf("%d, ",temp[i]);
    printf("\n");


    /* FINALIZE AND FREE */
    MPI_Finalize();
    if (rank==0)
        free(big_array);
    free(local_array);
    free(temp);
    return 0;
}
