#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include<unistd.h>

/*
Parallel quicksort using MPI
Assignment 3 - Parallel and Distributed Computing

By:
Claude Carlsson
Edvin Germer
Ture Hassler

2023-05-07
*/

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

// Merge lists function from internet
void merge(int* arr1, int* arr2, int arr1_size, int arr2_size, int* result) {
    int i = 0, j = 0, k = 0;
    while (i < arr1_size && j < arr2_size) {
        if (arr1[i] < arr2[j]) {
            result[k++] = arr1[i++];
        } else {
            result[k++] = arr2[j++];
        }
    }
    // Copy the remaining elements of arr1, if any
    while (i < arr1_size)
        result[k++] = arr1[i++];

    // Copy the remaining elements of arr2, if any
    while (j < arr2_size)
        result[k++] = arr2[j++];
}

// print list function
void print_list(int* list, int n)
{
    for (int i = 0; i < n; i++)
        printf("%d, ", list[i]);
    printf("\n\n");
}

int find_median(int *array, int n)
{
    int pivot;
    if(n%2==0)
        pivot = (array[(n-1)/2]+array[(n/2)])/2; // mean
    else
        pivot = array[n/2];

    return pivot;
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
            pivot = find_median(array, n);
        }
    } 
    else if (pivot_strategy == 2) // Median of medians
    {
        // find median
        pivot = find_median(array, n);
        
        // Allocate memory for all pivots
        if (rank == 0)
            pivots = (int*)malloc(size * sizeof(int));
        
        // gather pivot in pivots array on rank 0
        MPI_Gather(&pivot, 1, MPI_INT, pivots, 1, MPI_INT, 0, comm); // Gather all pivots to processor 0

        // sort pivots and find median of medians
        if (rank == 0)
        {
            quicksort(pivots,size);
            pivot = find_median(pivots, size);
        }
    } 
    else // Mean of medians  Pivot strategy 3
    {
        // find median
        pivot = find_median(array, n);
        
        // Allocate memory for all pivots
        if (rank == 0)
            pivots = (int*)malloc(size * sizeof(int));
        
        // gather pivot in pivots array on rank 0
        MPI_Gather(&pivot, 1, MPI_INT, pivots, 1, MPI_INT, 0, comm); // Gather all pivots to processor 0

        // find mean of medians
        if (rank == 0)
        {
            int sum=0;
            for (int i=0;i<size;i++)
                sum+=pivots[i];

            pivot = sum/size;
        }
    }

    // Broadcast pivot to all other processors in group
    MPI_Bcast(&pivot, 1, MPI_INT, 0, comm);

    // FREE MEMORY
    if (rank==0 && (pivot_strategy==2 || pivot_strategy==3))
        free(pivots);

    return pivot;
} 

// Inner loop of parallel quicksort algorithm
int par_quicksort(int **array_ptr, int n, int pivot_strategy, MPI_Comm comm)
{
    int pivot,rank,size,length;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Base case
    if (size < 2)
        return n;

    // Access list
    int *array = *array_ptr;
    
    // 3.1 Select Pivot element
    pivot = select_pivot(array, n, pivot_strategy, comm);

    // 3.2 Split the array into two subarrays and send to other processor 
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
    } 
    else // Send larger to the right:   
    {
        MPI_Recv(&recv_count, 1, MPI_INT, rank + size/2, 0, comm, MPI_STATUS_IGNORE);
        MPI_Send(&larger_count, 1, MPI_INT, rank + size/2, 0, comm);
    }

    // Send list elements
    int* temp = (int*)malloc(recv_count * sizeof(int));
    int* temp_keep; // part of array we are keeping (smaller or larger elements)

    if (rank >= size / 2) // Send smaller to the left:  0-2, 1-3
    {
        MPI_Send(array, smaller_count, MPI_INT, rank - size/2, 0, comm); 
        MPI_Recv(temp, recv_count, MPI_INT, rank - size/2, 0, comm, MPI_STATUS_IGNORE); 
        temp_keep = (int*)malloc(larger_count * sizeof(int));
        memcpy(temp_keep,&array[smaller_count], larger_count * sizeof(int)); // copy elements to keep

        // realloc for array
         length = recv_count + larger_count;
        if (length!=0)
        {
            array = (int*)realloc(array, (length) * sizeof(int));
        }
        // merge temp and temp_keep into array
        merge(temp_keep, temp, larger_count, recv_count, array);
  
    } 
    else // Send larger to the right:   
    {
        MPI_Recv(temp, recv_count, MPI_INT, rank + size/2, 0, comm, MPI_STATUS_IGNORE);
        MPI_Send(&array[smaller_count], larger_count, MPI_INT, rank + size/2, 0, comm);

        temp_keep = (int*)malloc(smaller_count * sizeof(int));
        memcpy(temp_keep, array, smaller_count * sizeof(int)); // copy elements to keep

        // realloc for array
        length = recv_count + smaller_count;
        if (length!=0)
        {
            array = (int*)realloc(array, (length) * sizeof(int));
        }
        // merge temp and temp_keep into array
        merge(temp_keep, temp, smaller_count, recv_count, array);
    }

    // free temp and temp_keep
    free(temp_keep);
    if (recv_count!=0)
        free(temp);

    // Split into groups
    MPI_Comm new_comm;
    MPI_Comm_split(comm, rank < size/2, rank, &new_comm);

    // Recursive quicksort call
    *array_ptr = array;
    int final_length = par_quicksort(array_ptr, length, pivot_strategy, new_comm);

    // Finalize and free
    MPI_Comm_free(&new_comm);
    return final_length;
}

// mpirun --bind-to none -n 4 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input125000000.txt output.txt 1
// mpirun --bind-to none -n 4 ./quicksort ./input16.txt output.txt 1
int main(int argc, char *argv[])
{
    /* INITIALIZE PARAMETERS */
    int rank,size,pivot,n, pivot_strategy;
    int *big_list, *local_list;

    /* INITIALIZE MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* CHECK ARGUMENT COUNT */
    if (rank==0)
    {
        if (argc != 4)
        {
            if (rank == 0)
                printf("Usage: %s <input_file> <output_file>\n", argv[0]);
            return -1;
        }
    }

    /* READ STRATEGY */
    pivot_strategy = atoi(argv[3]);

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
        // Close the input file
        fclose(input);
    }

    /* WAIT FOR DATA TO BE LOADED */
    MPI_Barrier(MPI_COMM_WORLD);

    /* DISTRIBUTE DATA */
    double start = MPI_Wtime(); // Record the start time
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // Broadcast the value of n to all processes

    // Determine local_list size
    int m = n / size;
    int remainder = n % size;
    if (rank < remainder)  // first ranks get remainder
        m = m + 1;
        
    int* send_counts, *displs; 
    if (rank == 0)
    {
        send_counts = (int*)malloc(size * sizeof(int)); // how many to send to each rank
        displs = (int*)malloc(size * sizeof(int));      // index where to start
        
        for (int i = 0; i < size; i++)
        {
            if (i < remainder)
                send_counts[i] = n/size+1;
            else
                send_counts[i] = n/size;

            if (i > 0)
                displs[i] = displs[i - 1] + send_counts[i - 1];
            else
                displs[i] = 0;
        }
    }

    // Allocate memory
    local_list = (int *)malloc(m* sizeof(int));

    // Scatter data
    MPI_Scatterv(big_list, send_counts, displs, MPI_INT, local_list, m, MPI_INT, 0, MPI_COMM_WORLD);

    /* SORT LOCALLY */
    quicksort(local_list, m);

    /* CALL PARALLEL QUICKSORT */
    int final_length = par_quicksort(&local_list, m, pivot_strategy, MPI_COMM_WORLD);

    /* GATHER ALL FINAL LENGTHS */
    int* lengths = NULL;
    if (rank == 0)
        lengths = (int*) malloc(size * sizeof(int));
    MPI_Gather(&final_length, 1, MPI_INT, lengths, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    /* GATHER SORTED LOCAL ARRAYS */
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

    /* LOCAL TIME */
	double local_execution_time = MPI_Wtime() - start;
    
    /* SAVE LIST TO OUTPUT FILE */
    /* if (rank == 0) {
        FILE* output = fopen(argv[2], "w");
        for (int i = 0; i < n; i++)
            fprintf(output, "%d ", big_list[i]);
            fclose(output);
    } */

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

    /* FINALIZE AND FREE MEMORY */
    if (rank==0)
    {
        free(big_list);
        free(lengths); 
        free(send_counts);
        free(displs);
    }
    if (final_length!=0)
        free(local_list);

    MPI_Finalize();
    return 0;
}