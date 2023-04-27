#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_matrix(double *matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.2lf ", matrix[i * cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// module load gcc openmpi
// mpirun --bind-to none -n 2 ./exercise2 1 3 4

// ./exercise2 iterations n m k
int main(int argc, char** argv)
{
    /* INITIALIZE PARAMETERS */
    int p, rank, iterations, n, m, k;

    int print=0;

    /* INITIALIZE MPI*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* CHECK INPUT */
    if (argc != 4)
    {
        if (rank == 0)
            printf("Incorrect number of input parameters");
        MPI_Finalize();
        return 1;
    }

    /* READ PARAMETERS */
    iterations = atoi(argv[1]);  // Times to repeat eq (2)
    n = atoi(argv[2]);           // Columns of matrix A
    m = atoi(argv[3]);           // Rows of matrix A
   // k = atoi(argv[4]);           // Size of blocks  (kxn)

    k = m/p; // automatically set k

    /* MEMORY ALLOCATION AND MATRIX GENERATION */
    double *A, *A_local, *Ax_local, *x, *w, *w_local;
    
    x = (double *)malloc(n * sizeof(double));
    A_local = (double *)malloc(k * n * sizeof(double));
    Ax_local = (double *)malloc(k * sizeof(double));
    w_local = (double *)malloc(n * sizeof(double));

    if (rank==0)
    {
        A = (double *)malloc(m * n * sizeof(double));
        w = (double *)malloc(n * sizeof(double));

        // Matrix A
        for (int i=0;i<m;i++)
            for (int j=0;j<n;j++)
                A[i*n + j] = rand()%10;
        
        // Vector X
        for (int i=0;i<n;i++)
            x[i] = rand()%10;

        if (print==1)
        {
            printf("Large A\n");
            print_matrix(A,m,n);
            printf("\nx\n");
            print_matrix(x,n,1);
            printf("____________________________\n");
        }
    }


    /* START TIMER */
    double start = MPI_Wtime();

    // DISTRIBUTE MATRIX AND VECTOR X
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Send x to all processes
    MPI_Scatter(A, k * n, MPI_DOUBLE, A_local, k * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); // send row-blocks of A to all processes

    if (print==1)
    {
        printf("Rank %d A:\n",rank);
        print_matrix(A_local,k,n);
        printf("_____________________\n");
    }
    
    /* CALCULATIONS */ 
    for (int iter=0;iter<iterations;iter++)
    {
        // Calculate Ax
        for (int i=0;i<k;i++)  // k rows in A_local
        {
            double res = 0.0;
            for (int j=0;j<n;j++)
            {
                res += A_local[i*n + j] * x[j];
            }
            Ax_local[i] = res;
        }

        // Calculate W
        for (int i=0;i<n;i++)  // n rows in A_local^T
        {
            double res = 0.0;
            for (int j=0;j<k;j++)
            {
                res += A_local[i + j*n] * Ax_local[j];
            }
            w_local[i] = res;
        }

        // Gather resuls from all processes
        MPI_Allreduce(w_local, x, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // Update x
        MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Send x to all processes
    }

    /* STOP TIMER */
    double local_execution_time = MPI_Wtime();

    /* FINAL W */
    if (print==1)
    {
        if (rank==0)
        {
            printf("__________________________________");
            printf("\n\nGathered final w:\n");
            print_matrix(x,n,1);

        }
    }

    /* FIND LONGEST TIME */
	double timings[p];
	MPI_Gather(&local_execution_time, 1, MPI_DOUBLE, timings, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		double max = 0;
		for (int i=0; i<p;i++)
		{
			double time = timings[i];
			if (time>max)
				max = time;
		}
		printf("%f\n", max);
	}

    /* FINALIZE AND FREE */
    if (rank==0)
    {
        free(A);
        free(w);
    }
    free(x);
    free(A_local);
    free(Ax_local);
    free(w_local);
   
    MPI_Finalize();
    return 0;
}
