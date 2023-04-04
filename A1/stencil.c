#include "stencil.h"

// Files:
// ls /home/maya/public/PDP_Assignment1/

/* 96 DATAPOINTS */
// mpirun --bind-to none -n 1 ./stencil /home/maya/public/PDP_Assignment1/input96.txt ./output.txt 1
// mpirun --bind-to none -n 2 ./stencil /home/maya/public/PDP_Assignment1/input96.txt ./output.txt 1
// mpirun --bind-to none -n 4 ./stencil /home/maya/public/PDP_Assignment1/input96.txt ./output.txt 1
// mpirun --bind-to none -n 8 ./stencil /home/maya/public/PDP_Assignment1/input96.txt ./output.txt 1

// mpirun --bind-to none -n 1 ./stencil /home/maya/public/PDP_Assignment1/input96.txt ./output.txt 4
// mpirun --bind-to none -n 2 ./stencil /home/maya/public/PDP_Assignment1/input96.txt ./output.txt 4
// mpirun --bind-to none -n 4 ./stencil /home/maya/public/PDP_Assignment1/input96.txt ./output.txt 4
// mpirun --bind-to none -n 8 ./stencil /home/maya/public/PDP_Assignment1/input96.txt ./output.txt 4


int main(int argc, char **argv)
{
	/* The program  reads the input file, applies the stencil the specified number of times
     and writes the result to the output file. */
	
	if (4 != argc) {
		printf("Usage: stencil input_file output_file number_of_applications\n");
		return 1;
	}
	char *input_name = argv[1];    // Input file name
	char *output_name = argv[2];   // Output file name
	int num_steps = atoi(argv[3]); // How many times the stencil will be applied

	int num_values;                // Total elements in input data
	double* input;
	double* output;
	int rank, size;

	/* SETUP MPI */
	MPI_Status status;

	MPI_Init(&argc, &argv);               /* Initialize MPI               */
	MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */


	/*  I/O HANDLED BY RANK 0 */
	if (rank==0)
	{
		// Read data
		if (0 > (num_values = read_input(input_name, &input)))  // TAKES A LONG TIME (40s for N=100.000.000)
		{
			return 2;
		}

		// Allocate data for result
		if (NULL == (output = malloc(num_values * sizeof(double)))) {
			perror("Couldn't allocate memory for output");
			return 2;
		}
	}


	/* SEND TOTAL NUMBER OF ELEMENTS TO ALL PROCESSES FROM RANK 0 */
	MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	
	/* DISTRIBUTE DATA FROM RANK 0 */
	int batch_size = num_values / size;
	double* local_input = malloc(batch_size * sizeof(double));
	double* local_output = malloc(batch_size * sizeof(double));
	MPI_Scatter(input, batch_size, MPI_DOUBLE, local_input, batch_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	/* CREATE STENCIL */
	double h = 2.0*PI/num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH/2;
	const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};


	/* START TIMER */
	double start = MPI_Wtime();
	

	/* APPLY STENCIL ON ALL PROCESSES*/
	double first[EXTENT];
	double last[EXTENT];
	int prev_rank = (rank - 1 + size) % size;
	int next_rank = (rank + 1) % size;

	MPI_Request send_reqs[2], recieve_reqs[2]; // 2 send and 2 recieve for each process

	for (int s=0; s<num_steps; s++) // Repeatedly apply stencil
	{
		if (size>1)
		{
			// Send to the previous rank
			MPI_Isend(&local_input[0], EXTENT, MPI_DOUBLE, prev_rank, rank, MPI_COMM_WORLD, &send_reqs[0]); // Send first 2
			MPI_Irecv(&last[0], EXTENT, MPI_DOUBLE, next_rank, next_rank, MPI_COMM_WORLD, &recieve_reqs[0]); // Recieve last 2
			
			// Send to the next rank
			MPI_Isend(&local_input[batch_size-EXTENT], EXTENT, MPI_DOUBLE, next_rank, rank, MPI_COMM_WORLD, &send_reqs[1]); // Send last 2
			MPI_Irecv(&first[0], EXTENT, MPI_DOUBLE, prev_rank, prev_rank, MPI_COMM_WORLD, &recieve_reqs[1]); // Recieve first 2
		}
		
		
		// Apply stencil on middle elements i
		for (int i=EXTENT; i<batch_size-EXTENT; i++) // Middle elements
		{
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) // Sweep stencil
			{
				int index = i - EXTENT + j;
				result += STENCIL[j] * local_input[index];
			}
			local_output[i] = result;
		}

		// Using multiple processes. Wait to recieve halo elements
		if (size>1)
			MPI_Waitall(2, recieve_reqs, MPI_STATUSES_IGNORE);


		// Apply stencil on first elements i
		for (int i=0; i<EXTENT; i++)  // First two elements
		{
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) // Sweep stencil 
			{
				if (size==1) // Default serial code
				{
					int index = (i - EXTENT + j + num_values) % num_values;
					result += STENCIL[j] * local_input[index];
				}
				else // Parallel code
				{
					int index = i - EXTENT + j;
					if (index < 0)
						result += STENCIL[j]*first[index+EXTENT];
					else
						result += STENCIL[j]*local_input[index];
				}
			}
			local_output[i] = result;
		}

		// Apply stencil on last elements i
		for (int i=batch_size-EXTENT; i<batch_size; i++) // Last two elements
		{
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) // Sweep stencil 
			{
				if (size==1) // Default serial code
				{
					int index = (i - EXTENT + j) % num_values;
					result += STENCIL[j] * local_input[index];
				}
				else // Parallel code
				{
					int index = i - EXTENT + j;
					if (index >= batch_size)
						result += STENCIL[j]*last[index-batch_size];
					else
						result += STENCIL[j]*local_input[index];
				}
			}
			local_output[i] = result;
		}
	

		// Swap input and output
		if (s < num_steps-1) {
			double *tmp = local_input;
			local_input = local_output;
			local_output = tmp;
		}
	}


	/* LOCAL TIME */
	double local_execution_time = MPI_Wtime() - start;


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
		printf("      LONGEST TOTAL TIME: %f\n", max);
	}


	/* FINAL COLLECTION - SEND DATA BACK TO RANK 0 */
	MPI_Gather(local_output, batch_size, MPI_DOUBLE, output, batch_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	
	/*  I/O HANDLED BY RANK 0 */
	if (rank == 0)
	{ 
		#ifdef PRODUCE_OUTPUT_FILE
		if (0 != write_output(output_name, output, num_values)) // Write result
			return 2;
		#endif
	}
	 
	if(rank==0)
	{
		free(output);
		free(input);
	}
	free(local_input);
	free(local_output);

	/* CLEAN UP*/
	MPI_Finalize();

	return 0;
}


/* FUNCTION FOR READING INPUT FILE*/
/*          DO NOT CHANGE         */
int read_input(const char *file_name, double **values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values)) {
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num_values;
}

/* FUNCTION FOR CREATING OUTPUT FILE*/
/*          DO NOT CHANGE           */
int write_output(char *file_name, const double *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%.4f ", output[i])) {
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
	}
	return 0;
}
 