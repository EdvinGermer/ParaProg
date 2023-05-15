#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <unistd.h>

/*
Monte Carlo Computatio for Malaria Epidemic

By:
Edvin Germer
2023-05-12
*/

// Header for prop function
void prop(int *x, double *w);



// mpirun --bind-to none -n 4 ./montecarlo 16
int main(int argc, char *argv[])
{   
    /* INITIALIZE PARAMETERS */
    int rank,size;
    int N,n;
    double t;     // Discrete time variable
    double tau;   // Time step
    int T=100; // Total simulation time per experiment
    int b=20;  // Number of histogram bins

    double a0;
    double u1;   // Random number 1
    double u2;   // Random number 2

    double w[15];


    /* INITIALIZE MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* READ INPUTS */
    if (argc != 2)
    {
        if (rank == 0)
            printf("Usage:  ./montecarlo N\n");
        return -1;
    }
    
    /* HOW MANY EXPERIMENTS PER PROCESS */
    N = atoi(argv[1]); // Total number of experiments
    n = N/size;           // Experiments per process

    /* DEFINE MATRIX P */
    int P[15][7] = {
    {1, 0, 0, 0, 0, 0, 0},
    {-1, 0, 0, 0, 0, 0, 0},
    {-1, 0, 1, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0, 0},
    {0, -1, 0, 0, 0, 0, 0},
    {0, -1, 0, 1, 0, 0, 0},
    {0, 0, -1, 0, 0, 0, 0},
    {0, 0, -1, 0, 1, 0, 0},
    {0, 0, 0, -1, 0, 0, 0},
    {0, 0, 0, -1, 0, 1, 0},
    {0, 0, 0, 0, -1, 0, 0},
    {0, 0, 0, 0, -1, 0, 1},
    {0, 0, 0, 0, 0, -1, 0},
    {1, 0, 0, 0, 0, 0, -1},
    {0, 0, 0, 0, 0, 0, -1}
    };

    /* START LOCAL TIME */
    double start = MPI_Wtime();

    /* MONTE CARLO SIMULATION */
    int **X = calloc(7, sizeof(int*));
    for (int i=0; i<7; i++)
        X[i] = calloc(n, sizeof(int));
    
    
    srand(time(NULL)+rank); // Random seed for each rank

    for (int i=0;i<n;i++)  
    {
        int x[7] = {900,900,30,330,50,270,20};  // Initial state vector
        
        // Do experiment (simulation)
        t=0;
        while (t<T)
        {
            // Compute w
            prop(x, w);
            
            // Compute a0
            a0 = 0;
            for (int i=0;i<15;i++)
                a0 += w[i];

            // Generate two random numbers
            u1 = (double)rand() / RAND_MAX;
            u2 = (double)rand() / RAND_MAX;

            // Set time increment tau
            tau = -log(u1)/a0;

            // Find r
            double sum=0;
            double target = a0*u2;
            int r;
            
            for (r=0;r<15;r++)
            {
                sum += w[r];
                if (sum>=target)
                    break;
            }

            // Update state vector 
            for (int i=0;i<7;i++)
                x[i] += P[r][i];

            // Update time
            t += tau;
        }

        // Save results
        for (int row=0;row<7;row++)
            X[row][i] = x[row];   // save x as a column in X
    }

    /* SANITY CHECK */
    /* printf("Rank %d: \n", rank);
    for (int i=0;i<n;i++)
        printf("%d, ", X[0][i]);
    printf("\n"); */


    /* EXTRACT SUSCEPTIBLE HUMANS */
    int* local_res = (int*)malloc(n * sizeof(int));

    for (int i=0;i<n;i++)
        local_res[i] =  X[0][i];


    /* FIND LOCAL MAX AND MIN */
    int local_min = local_res[0];
    int local_max = local_res[0];

    for (int i=1;i<n;i++)
    {
        if (local_res[i]<=local_min)
            local_min = local_res[i];
        if (local_res[i]>=local_max)
            local_max = local_res[i];
    }
    //printf("RANK %d:   local_max = %d, local_min = %d\n", rank,local_max,local_min);


    /* FIND GLOBAL MAX AND MIN */
    int global_max, global_min;

    MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    

    /* SANITY CHECK */
    //printf("RANK %d:  local_max = %d, local_min = %d\n",rank, local_max, local_min);
    //if (rank==0)
    //    printf("RANK %d:  global_max = %d, global_min = %d\n",rank, global_max, global_min);


    /* LOCAL HISTOGRAM */
    double interval = (global_max-global_min)/b;  // width of each bar

    int bins[b+1];  // array of bar edges
    for (int i=0;i<b;i++)
        bins[i] = global_min + i*interval;
    bins[b] = global_max;
    

    int freqs[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // count each bin
    int bin;
    for (int i=0;i<n;i++)
    {
        bin = (local_res[i]-global_min)/interval; // determine which bin
        freqs[bin]++; // increment that bin
    }


    /* GATHER HISTOGRAM */
    int histogram[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    MPI_Reduce(freqs, histogram, b, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


    /* STOP LOCAL TIME */
	double local_execution_time = MPI_Wtime() - start;


    /* WRITE RESULTS TO OUTPUT FILE */
    if (rank == 0)
    {
        FILE *file = fopen("histogram.txt", "w");

        // Print histogram bins
        for (int i = 0; i <= b; i++)
        {
            fprintf(file, "%d", bins[i]);
            if (i < b)
                fprintf(file, ", ");
        }

        // Newline
        fprintf(file, "\n");

        // Print histogram freqs
        for (int i = 0; i < b; i++)
        {
            fprintf(file, "%d", histogram[i]);
            if (i < b - 1)
                fprintf(file, ", ");
        }
        fclose(file);
    }


    /* FIND LONGEST TIME */
    double max_time;
    MPI_Reduce(&local_execution_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (rank == 0)
	    printf("%f\n", max_time);


    /* FINALIZE AND FREE MEMORY */
    free(local_res);
    
    for (int i=0; i<7; i++)
        free(X[i]);
    free(X);

    MPI_Finalize();
    return 0;
}