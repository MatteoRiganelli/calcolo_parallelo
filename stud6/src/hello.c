#include <mpi.h>
#include <stdio.h>

int main(argc, argv)
	int argc;
	char *argv[]; {
	int size, rank, rc;
	
	rc=MPI_Init(&argc, &argv);
	if(rc!=MPI_SUCCESS) {
		printf("Inizialization error.\n");
		MPI_ABORT(MPI_COMM_WORLD, &size);
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	printf("Num di task=%d - rank=%d\n", size, rank);

	MPI_Finalize();

}
