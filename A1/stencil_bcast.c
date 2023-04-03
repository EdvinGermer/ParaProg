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
	

	/* SETUP MPI */
	MPI_Status status;
	int rank, size;

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


	/* SEND TOTAL NUMBER OF ELEMENTS TO ALL PROCESSES FROMT RANK 0 */
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
	double stencil_start, stencil_end, com_end, com_start;

	double communication_time = 0.0;
	double stencil_time = 0.0;


	/* BROADCAST HALO ELEMENTS */
	double halo[4*size]; // Array for storing halo elements
	if (rank==0)
	{
		for(int i=0;i<EXTENT*2*size;i+=4)
		{
			int idx = i*batch_size;

			halo[i] = input[idx];
			halo[i + 1] = input[idx+1];

			halo[i + 2] = input[idx+batch_size-2];
			halo[i + 3] = input[idx+batch_size+-1];
		}
	}
	com_start = MPI_Wtime();
	MPI_Bcast(&halo, 4*size, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Send to all
	com_end = MPI_Wtime();
	communication_time += com_end-com_start;


	/* FIND EACH RANKS HALO ELEMENTS */
	double first[EXTENT];
	double last[EXTENT];

	int first_halo_idx = (rank * 4 + 4 * size - 2) % (4 * size);
	int last_halo_idx = (rank * 4 + 2) % (4 * size);

	first[0] = halo[first_halo_idx];
	first[1] = halo[first_halo_idx + 1];

	last[0] = halo[last_halo_idx];
	last[1] = halo[last_halo_idx + 1];


	/* APPLY STENCIL ON ALL PROCESSES*/
	stencil_start = MPI_Wtime();
	for (int s=0; s<num_steps; s++) // Repeatedly apply stencil
	{
		// Apply stencil on elements i
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
	stencil_end = MPI_Wtime();
	stencil_time += stencil_end-stencil_start;

	/* LOCAL TIME */
	double local_execution_time = MPI_Wtime() - start;

	//printf("    Rank %d took %f seconds\n",rank, local_execution_time);
	//printf("    	Com_time = %f, stencil_time = %f\n",communication_time, stencil_end-stencil_time);


	/* FIND LONGEST TIME */
	double timings0[size];
	double timings1[size];
	double timings2[size];
	MPI_Gather(&local_execution_time, 1, MPI_DOUBLE, timings0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&communication_time, 1, MPI_DOUBLE, timings1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(&stencil_time, 1, MPI_DOUBLE, timings2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		double max0 = 0;
		double max1 = 0;
		double max2 = 0;
		for (int i=0; i<size;i++)
		{
			double time0 = timings0[i];
			double time1 = timings1[i];
			double time2 = timings2[i];
			if (time0>max0)
				max0 = time0;
			if (time1>max1)
				max1 = time1;
			if (time2>max2)
				max2 = time2;
		}
		//printf("    LONGEST TOTAL TIME = %f\n", max0);
		//printf("    LONGEST COMMUNICATION TIME = %f\n", max1);
		//printf("    LONGEST STENCIL TIME = %f\n", max2);
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

	
	/* CLEAN UP*/
	MPI_Finalize();
	
	if(rank==0)
	{
		free(output);
		free(input);
	}
	free(local_input);
	free(local_output);

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
