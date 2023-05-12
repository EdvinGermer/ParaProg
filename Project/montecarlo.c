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



// mpirun --bind-to none -n 4 ./montecarlo 1000
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
    double* local_res;


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

    /* DEFINE VECTOR X AND MATRIX P */
    int x[7] = {900,900,30,330,50,270,20};  // Initial state vector

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
    local_res = (double *)malloc(n * sizeof(double));
    srand(time(NULL)); // Random seed

    for (int i=0;i<n;i++)  
    {
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
            double sum1=0,sum2=0;
            int r = 0;
            double target = a0*u2;
            int counter;
            while (counter<15)
            {
                // sum1
                for (int i=0;i<(r-1);i++)
                    sum1 += w[i];
                // sum2
                for (int i=0;i<r;i++)
                    sum2+= w[i];

                // check condition
                if (sum1<target && sum2>=target) // exit loop with current value on r
                    counter=15; 

                // increment
                counter++;
                r++;
            }

            // Update state vector 
            for (int i;i<7;i++)
                x[i] += P[r][i];

            // Update time
            t += tau;
        }

        // Save results
        local_res[i] = a0;
    }

    /* GATHER RESULTS */


    /* STOP LOCAL TIME */
	double local_execution_time = MPI_Wtime() - start;

    /* FIND LONGEST TIME */
    double max;
    MPI_Reduce(&local_execution_time, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (rank == 0)
	    printf("%f\n", max);

    /* FINALIZE AND FREE MEMORY */
    MPI_Finalize();
    return 0;
}