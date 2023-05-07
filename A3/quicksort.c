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

int select_pivot(int *array, int n, int pivot_strategy, MPI_Comm comm)
{
    int size,rank,pivot;
    int* pivots;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    if (pivot_strategy == 1) // Median of processor 0 
    {
        if (rank == 0) // Select the median on rank 0
        {
            if(n%2==0)
                pivot = (array[(n-1)/2]+array[(n/2)])/2; // mean
            else
                pivot = array[n/2];
        }
    } 
    else if (pivot_strategy == 2) // Median of medians
    {
        pivot = array[n/2];
        MPI_Gather(&pivot, 1, MPI_INT, &pivot, 1, MPI_INT, 0, comm); // Gather all pivots to processor 0

        if (rank == 0)
            pivot = array[n/2]; // Select the median

    } else // Mean of medians
    {
        pivot = array[n/2];

        if (rank == 0)
            pivots = (int*)malloc(size * sizeof(int)); // Allocate memory for all pivots

        MPI_Gather(&pivot, 1, MPI_INT, pivots, size, MPI_INT, 0, comm); // Gather all pivots to processor 0
        
        if (rank == 0) // Mean of all medians 
        { 
            int sum = 0;
            for (int i = 0; i < size; i++)
                sum += pivots[i];
            pivot = sum/size;
        }
    }

    // Broadcast pivot to all other processors in group
    MPI_Bcast(&pivot, 1, MPI_INT, 0, comm);
    //printf("Rank %d, Pivot: %d\n",rank, pivot);

    // FREE MEMORY
    if (rank==0 && (pivot_strategy==2 || pivot_strategy==3))
        free(pivots);

    return pivot;
} 

int par_quicksort(int **array_ptr, int n, int pivot_strategy, MPI_Comm comm)
{
    int pivot;
    int rank;
    int size;
    int length;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Base case
    if (size < 2)
        return n;

    // Access list
    int *array = *array_ptr;

    //if (rank==0)
        //printf("------- GROUPING -------\n");
    
    // 3.1 Select Pivot element
    pivot = select_pivot(array, n, pivot_strategy, comm);
    //if (rank==0)
        //printf("Pivot element = %d\n",pivot);

    /// 3.2 Split the array into two subarrays and send to other processor 
    int smaller_count = 0;
    for (int i = 0; i < n; i++)
        if (array[i] <= pivot)
            smaller_count++;

    int larger_count = n-smaller_count;

    // Exchange send_count
    int recv_count = 0;
    if (rank >= size / 2) // Send smaller to the left:  0-2, 1-3
    {
        MPI_Send(&smaller_count, 1, MPI_INT, rank - size/2, 0, comm); // Send nr of smaller elements
        MPI_Recv(&recv_count, 1, MPI_INT, rank - size/2, 0, comm, MPI_STATUS_IGNORE); // Receive recv_count
        //printf("Rank %d with pair %d: recieved %d and sent %d\n",rank,rank - size/2,recv_count,smaller_count);
    } 
    else // Send larger to the right:   
    {
        MPI_Recv(&recv_count, 1, MPI_INT, rank + size/2, 0, comm, MPI_STATUS_IGNORE);
        MPI_Send(&larger_count, 1, MPI_INT, rank + size/2, 0, comm);
        //printf("Rank %d with pair %d: recieved %d and sent %d\n",rank,rank + size/2,recv_count,larger_count);
    }

    // Send list elements
    int* temp = (int*)malloc(recv_count * sizeof(int));
    if (rank >= size / 2) // Send smaller to the left:  0-2, 1-3
    {
        MPI_Ssend(array, smaller_count, MPI_INT, rank - size/2, 0, comm); 
        MPI_Recv(temp, recv_count, MPI_INT, rank - size/2, 0, comm, MPI_STATUS_IGNORE); 
        memcpy(&array[0],&array[smaller_count], larger_count * sizeof(int)); // move elements to start
        length = larger_count + recv_count;
        if (length!=0)
            {
                array = (int*)realloc(array, (length) * sizeof(int));
                memcpy(&array[larger_count],temp, recv_count * sizeof(int)); // append temp to end
            }
        //printf("Rank %d:\n",rank);
        //print_list(array,larger_count + recv_count);
    } 
    else // Send larger to the right:   
    {
        MPI_Recv(temp, recv_count, MPI_INT, rank + size/2, 0, comm, MPI_STATUS_IGNORE);
        MPI_Ssend(&array[smaller_count], larger_count, MPI_INT, rank + size/2, 0, comm);
        length = smaller_count + recv_count;
        if (length!=0)
            {
                array = (int*)realloc(array, (length) * sizeof(int)); // realloc at end of array
                memcpy(&array[smaller_count],temp, recv_count * sizeof(int)); // append temp to end
            }
        //printf("Rank %d:\n",rank);
        //print_list(array,smaller_count + recv_count);
    }



    // Local Sort
    quicksort(array, length);
    //printf("Rank %d:\n",rank);
    //print_list(array,length);

    //if (length==0)
    //    printf("rank %d\n",rank);

    // Split into groups
    MPI_Comm new_comm;
    MPI_Comm_split(comm, rank < size/2, rank, &new_comm);

    // Recursive quicksort call
    *array_ptr = array;
    int final_length = par_quicksort(array_ptr, length, pivot_strategy, new_comm);

    // Finalize and free
    MPI_Comm_free(&new_comm);
    if (recv_count!=0)
        free(temp);
    return final_length;
}


// mpirun --bind-to none -n 4 ./quicksort ./input16.txt output.txt

int main(int argc, char *argv[])
{
    /* INITIALIZE PARAMETERS */
    int rank,size,pivot;
    int n;

    int* big_list;
    int* local_list;

    /* INITIALIZE MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the value of n to all processes

    // Determine local_list size
    int m = n/size;

    // Allocate memory
    local_list = (int *)malloc(m* sizeof(int));

    // Scatter data
    MPI_Scatter(big_list, m, MPI_INT, local_list, m, MPI_INT, 0, MPI_COMM_WORLD);

    // Print local_list before sort
    /* printf("RANK %d: ", rank);
    print_list(local_list,m); */




    /* SORT LOCALLY */
    quicksort(local_list,m);

    // Print local_list after local sort
    /* printf("RANK %d sorted: ", rank);
    print_list(local_list,m); */





    /* CALL PARALLEL QUICKSORT */
    int pivot_strategy = 1;
    int final_length = par_quicksort(&local_list, m, pivot_strategy, MPI_COMM_WORLD);




    /* GATHER ALL FINAL LENGTHS */
    int* lengths = NULL;
    if (rank == 0)
        lengths = (int*) malloc(size * sizeof(int));
    MPI_Gather(&final_length, 1, MPI_INT, lengths, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //printf("RANK %d, final_length = %d\n",rank, final_length);
    //print_list(local_list,final_length);
    

    /* PRINT LOCAL LISTS */
    printf("_________________________________\nRANK %d: SORTED LOCAL_LIST:\n",rank);
    print_list(local_list, final_length);



    /* GATHER SORTED LOCAL ARRAYS */
    int *displs;
    if (rank == 0)
    {
        // Find displacements
        displs = (int *)malloc(size * sizeof(int));
        displs[0] = 0;
        for (int i = 1; i < size; i++)
            displs[i] = displs[i - 1] + lengths[i - 1];
    
        // Save local_lists to big_list
        memcpy(big_list, local_list, final_length * sizeof(int)); // Copy rank 0 data directly
        for (int i = 1; i < size; i++)
            MPI_Recv(big_list + displs[i], lengths[i], MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else // Send local_list to rank 0
        MPI_Send(local_list, final_length, MPI_INT, 0, 0, MPI_COMM_WORLD);


    /* PRINT FINAL SORTED LIST*/
    if (rank == 0)
    {
        printf("_________________________________\nSORTED BIG LIST:\n");
        print_list(big_list, n);
    }

    /* SAVE LIST TO OUTPUT FILE */
    if (rank == 0) {
        FILE* output = fopen(argv[2], "w");
        for (int i = 0; i < n; i++)
            fprintf(output, "%d ", big_list[i]);
            fclose(output);
    }

    /* FINALIZE AND FREE MEMORY */
    if (rank==0)
    {
        free(big_list);
        free(lengths);   
        free(displs);
    }
    //if (final_length!=0)
    //    free(local_list);

    MPI_Finalize();
    return 0;
}