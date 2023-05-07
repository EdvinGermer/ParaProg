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

void par_quicksort(int *array, int n, int pivot_strategy, MPI_Comm comm)
{
    printf("Started par_quicksort\n");

    int pivot;
    int rank;
    int size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Oklart om det här behövs
    if (size < 2)
        return;

    // 3.1 Select Pivot element ()
    if (pivot_strategy == 1) // Median of processor 0 
    {
        if (rank == 0)
        {
            // Select the median
            if(n%2==0)
            {
                pivot = (array[n/2]+array[(n/2)+1])/2; // mean
            }
            else
            {
                pivot = array[n/2]
            }
        }
    } else if (pivot_strategy == 2) // Median of medians
    {
        pivot = array[n/2];
        MPI_Gather(&pivot, 1, MPI_INT, &pivot, 1, MPI_INT, 0, comm); // Gather all pivots to processor 0

        if (rank == 0)
        {
            // Select the median
            pivot = array[n/2];
        }


    } else // Mean of medians
    {
        pivot = array[n/2];
        int* pivots;

        if (rank == 0)
        {
            // Allocate memory for all pivots
            pivots = (int*)malloc(size * sizeof(int));
        }

        MPI_Gather(&pivot, 1, MPI_INT, &pivots, size, MPI_INT, 0, comm); // Gather all pivots to processor 0
        
        if (rank == 0)
        {
            // Mean of all medians 
            int sum = 0;
            for (int i = 0; i < size; i++)
                sum += pivots[i];
            pivot = sum/size;
        }
    }

    printf("After pivot selection\n");
    printf("Pivot: %d\n", pivot);

    // Broadcast pivot to all other processors in group
    MPI_Bcast(&pivot, 1, MPI_INT, 0, comm);

    // 3.2 Split the array into two subarrays and send to other processor 
    int* smaller = (int*)malloc(n * sizeof(int));
    int* larger = (int*)malloc(n * sizeof(int));
    int smaller_size = 0;
    int larger_size = 0;

    for (int i = 0; i < size; i++)
    {
        if (array[i] < pivot)
        {
            smaller[smaller_size] = array[i];
            smaller_size++;
        } else
        {
            larger[larger_size] = array[i];
            larger_size++;
        }
    }

    printf("After split\n");
    printf("Now sending to other processor\n");

    if (rank > size / 2) // Send larger to other processor, 1-8, 2-7, 3-6, 4-5
    {
        // Send size of small and then the array
        MPI_Isend(smaller_size, 1, MPI_INT, rank - size/2, 0, comm);
        MPI_Isend(smaller, smaller_size, MPI_INT, rank - size/2, 0, comm);

        // Receive size of large and then the array
        MPI_Irecv(larger_size, 1, MPI_INT, rank - size/2, 0, comm, MPI_STATUS_IGNORE);
        MPI_Irecv(larger, smaller_size, MPI_INT, rank - size/2, 0, comm, MPI_STATUS_IGNORE);
    } else // Send smaller to other processor, 1-2, 2-3, 3-4, 4-5
    {
        MPI_Isend(larger_size, 1, MPI_INT, rank + size/2, 0, comm);
        MPI_Irecv(smaller, smaller_size, MPI_INT, rank + size/2, 0, comm, MPI_STATUS_IGNORE);
        
    }

    printf("After send and receive\n");
    printf("Splitting into two groups\n");

    // Split into two MPI groups
    MPI_Comm smaller_comm;
    MPI_Comm larger_comm;
    MPI_Comm_split(comm, rank < size/2, rank, &smaller_comm);
    MPI_Comm_split(comm, rank >= size/2, rank, &larger_comm);

    // 3.3 Recursively call par_quicksort on the two subarrays
    par_quicksort(smaller, smaller_size, pivot_strategy, smaller_comm);
    par_quicksort(larger,  larger_size, pivot_strategy, larger_comm);

    // Oklart hur man slår ihop listorna?? 
    // Gör det här? MPI Gather? If rank==0??
    
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

    // Do parallel quicksort 
    printf("Starting parallel sort: \n");
    par_quicksort(local_list, m, 1, MPI_COMM_WORLD);

    // Print local_list
    printf("RANK %d parallel sorted: ", rank);
    print_list(local_list,m);



   

    MPI_Finalize();
    return 0;
}

