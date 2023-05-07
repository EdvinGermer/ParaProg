#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

// Function to swap two numbers
void swap(double *arr, int i, int j)
{
    double t = arr[i];
    arr[i] = arr[j];
    arr[j] = t;
}

// Function to calculate the median of an array
double calculate_median(double *arr, int start, int end)
{
    int n = end - start + 1;
    int median;

    if (n % 2 == 0)
    {
        median = (arr[start + (n / 2) - 1] + arr[start + (n / 2)]) / 2;
    }
    else
    {
        median = arr[start + (n / 2)];
    }

    return median;
}

// Function to calculate the mean of an array
double calculate_mean(double *arr, int start, int end)
{
    int n = end - start + 1;
    int sum = 0;

    for (int i = start; i <= end; i++)
    {
        sum += arr[i];
    }

    return sum / n;
}

// Function that performs the Quick Sort
// for an array arr[] starting from the
// index start and ending at index end
// Function that performs the Quick Sort
// for an array arr[] starting from the
// index start and ending at index end
void quicksort(double *arr, int start, int end, int pivot_strategy, int rank, int size)
{
    if (start >= end)
        return;

    int pivot_index = start;
    double pivot;

    if (pivot_strategy == 1)
    {
        pivot_index = start + (end - start) / 2;
    }
    else if (pivot_strategy == 2 || pivot_strategy == 3)
    {
        double local_pivot;
        if (pivot_strategy == 2)
        {
            local_pivot = calculate_median(arr, start, end);
        }
        else 
        {
            local_pivot = calculate_mean(arr, start, end);
        }

        double *all_pivots = NULL;
        if (rank == 0)
        {
            all_pivots = (double *)malloc(size * sizeof(double));
        }

        MPI_Gather(&local_pivot, 1, MPI_DOUBLE, all_pivots, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            quicksort(all_pivots, 0, size - 1, 1, rank, size);
            pivot = all_pivots[size / 2];
            free(all_pivots);
        }

        MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int i = start; i <= end; i++)
        {
            if (arr[i] == pivot)
            {
                pivot_index = i;
                break;
            }
        }
    }

    pivot = arr[pivot_index];
    swap(arr, start, pivot_index);

    int i = start + 1;
    int j = start + 1;

    for (; j <= end; j++)
    {
        if (arr[j] < pivot)
        {
            swap(arr, i, j);
            i++;
        }
    }

    swap(arr, start, i - 1);

    quicksort(arr, start, i - 2, pivot_strategy, rank, size);
    quicksort(arr, i, end, pivot_strategy, rank, size);
}


// Function that merges the two arrays
double *merge(double *arr1, int n1, double *arr2, int n2)
{
    double *result = (double *)malloc((n1 + n2) * sizeof(double));
    int i = 0;
    int j = 0;

    for (int k = 0; k < n1 + n2; k++)
    {
        if (i >= n1)
        {
            result[k] = arr2[j];
            j++;
        }
        else if (j >= n2)
        {
            result[k] = arr1[i];
            i++;
        }
        else if (arr1[i] < arr2[j])
        {
            result[k] = arr1[i];
            i++;
        }
        else
        {
            result[k] = arr2[j];
            j++;
        }
    }
    return result;
}


// Driver Code
int main(int argc, char *argv[])
{
    int number_of_elements;
    double *data = NULL;
    int chunk_size, own_chunk_size;
    double *chunk;
    FILE *file = NULL;
    MPI_Status status;
    int rank, size;
    double *other_chunk;
    int other_chunk_size;
    int pivot_strategy;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 4)
    {
        if (rank == 0)
        {
            printf("Usage: mpirun -np <number_of_processors> %s <input_file> <output_file> <pivot_strategy>\n", argv[0]);
            printf("Pivot strategies:\n1 - Median in one processor\n2 - Median of all medians in each processor group\n3 - Mean value of all medians in each processor group\n");
        }
        MPI_Finalize();
        exit(0);
    }

    pivot_strategy = atoi(argv[3]);

    if (rank == 0)
    {
        file = fopen(argv[1], "r");
        if (file == NULL)
        {
            perror("File not found");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fscanf(file, "%d", &number_of_elements);
        data = (double *)malloc(number_of_elements * sizeof(double));

        for (int i = 0; i < number_of_elements; i++)
            fscanf(file, "%lf", &data[i]);

        fclose(file);

        chunk_size = number_of_elements / size;
    }

    MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&number_of_elements, 1, MPI_INT, 0, MPI_COMM_WORLD);

    own_chunk_size = chunk_size;
    if (rank == size - 1)
    {
        own_chunk_size = number_of_elements - chunk_size * (size - 1);
    }
    chunk = (double *)malloc(own_chunk_size * sizeof(double));

    MPI_Scatter(data, chunk_size, MPI_DOUBLE, chunk, chunk_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    quicksort(chunk, 0, own_chunk_size - 1, pivot_strategy, rank, size);

    for (int step = 1; step < size; step = 2 * step)
    {
        if (rank % (2 * step) != 0)
        {
            MPI_Send(&own_chunk_size, 1, MPI_INT, rank - step, 0, MPI_COMM_WORLD);
            MPI_Send(chunk, own_chunk_size, MPI_DOUBLE, rank - step, 0, MPI_COMM_WORLD);
            break;
        }

        if (rank + step < size)
        {
            MPI_Recv(&other_chunk_size, 1, MPI_INT, rank + step, 0, MPI_COMM_WORLD, &status);
            other_chunk = (double *)malloc(other_chunk_size * sizeof(double));
            MPI_Recv(other_chunk, other_chunk_size, MPI_DOUBLE, rank + step, 0, MPI_COMM_WORLD, &status);

            double *merged_array = merge(chunk, own_chunk_size, other_chunk, other_chunk_size);

            free(chunk);
            free(other_chunk);
            chunk = merged_array;
            own_chunk_size += other_chunk_size;
        }
    }

    if (rank == 0)
    {
        file = fopen(argv[2], "w");
        for (int i = 0; i < number_of_elements; i++) {
            double num = chunk[i];
            if (floor(num) == num && num != 0.0) {
                fprintf(file, "%.0lf\n", num);
            } else {
                fprintf(file, "%lf\n", num);
            }
        }
        fclose(file);
    }


    free(chunk);
    if (rank == 0)
        free(data);

    MPI_Finalize();
    return 0;
}
