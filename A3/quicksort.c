#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include<unistd.h>

/* Local quicksort function from the internet */ 
void quicksort(int *list, int n) 
{
    // Initialize
    int i, j, pivot, temp;

    // Base case
    if (n < 2)
        return;

    // Else
    pivot = list[n / 2]; 

    for (i = 0, j = n - 1;; i++, j--)
    {
        while (list[i] < pivot)
            i++;
        while (pivot < list[j])
            j--;
        if (i >= j)
            break;

        temp = list[i];
        list[i] = list[j];
        list[j] = temp;
    }
    quicksort(list, i);
    quicksort(list + i, n - i);
}

// print list function
void print_list(int* list, int n)
{
    for (int i = 0; i < n; i++)
        printf("%d, ", list[i]);
    printf("\n\n");
}


// mpirun --bind-to none -n 2 ./quicksort ./input4.txt output.txt

int main(int argc, char *argv[])
{
    /* INITIALIZE PARAMETERS */
    int rank,pivot;
    int n;

    int* big_list;
    int* local_list;

    /* INITIALIZE MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &pivot);

     /* CHECK ARGUMENT COUNT */
    if (rank==0)
    {
        if (argc != 3)
        {
            if (rank == 0)
                printf("Usage: %s <input_file> <output_file>\n", argv[0]);
            return -1;
        }
    }

    /* READ DATA TO RANK 0 */
    if (rank == 0)
    {
        FILE *input = fopen(argv[1], "r");

        // Read the size of the list
        if (fscanf(input, "%d", &n) != 1)
        {
            printf("Error when reading n from input file");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Allocate memory for list
        big_list = (int *)malloc(n * sizeof(int));

        // Read input into list
        for (int i = 0; i < n; i++)
        {
            if (fscanf(input, "%d", &big_list[i]) != 1)
            {
                printf("Error when reading A from input file");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // Print big_list
        printf("BIG LIST\n");
        print_list(big_list,n);

        // Close the input file
        fclose(input);
    }


    /* DISTRIBUTE DATA */
    // Broadcast the value of n to all processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Determine local_list size
    int m = n/pivot;

    // Allocate memory
    local_list = (int *)malloc(m* sizeof(int));

    // Scatter data
    MPI_Scatter(big_list, m, MPI_INT, local_list, m, MPI_INT, 0, MPI_COMM_WORLD);

    // Print local_list
    printf("RANK %d: ", rank);
    print_list(local_list,m);
 



    /* SORT LOCALLY */
    quicksort(local_list,m);

    // Print local_list
    printf("RANK %d sorted: ", rank);
    print_list(local_list,m);

   

    MPI_Finalize();
    return 0;
}