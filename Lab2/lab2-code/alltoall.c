#include        <stdio.h>
#include        <mpi.h>

void main(int argc, char *argv[]) {

	int ierror, rank, size;

	int nrecvd, sum, i, index, flag;
	MPI_Status send_status[10], recv_status[10];
	MPI_Request request[10];
	int out[10], in[10];

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (size > 10)
		MPI_Abort(MPI_COMM_WORLD, -1);

	for (i = 0; i < size; i++) {
		out[i] = rank + i;
		MPI_Isend(&out[i], 1, MPI_INT, i, 7, MPI_COMM_WORLD, &request[i]);
	}

	for (i = 0; i < size; i++)
		MPI_Recv(&in[i], 1, MPI_INT, i, 7, MPI_COMM_WORLD, &recv_status[i]);

	nrecvd = 0;
	while (nrecvd < size) {
		flag = 0;
		MPI_Testany(size, request, &index, &flag, send_status);

		if (flag) {
			nrecvd++;
			printf("%d finished sending to %d\n", rank, index);
		}
		else
			printf("%d waiting ...\n", rank);

	}

	printf("PE%d: in = ", rank);
	for (i = 0; i < size; i++)
		printf("%d ", in[i]);
	printf("\n");

	MPI_Finalize();
}
