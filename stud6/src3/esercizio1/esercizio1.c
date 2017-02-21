#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(argc, argv)

    int argc;
    char *argv[];   {

	double w_time = MPI_Wtime();

    int size, rank, rc,dest, somma=0;
    long int prodotto=1;

    MPI_Status status; /*Variabile di stato*/


    rc=MPI_Init(&argc, &argv);

    if(rc!=MPI_SUCCESS) {
           	printf("Inizialization error.\n");
            MPI_ABORT(MPI_COMM_WORLD, &size);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

   	int A;

    /* Intializes random number generator */
    srand(time(NULL) + rank);

    //for(j=1;j<num;j++) {
    A=(rand() % 100) +1;
    printf("my rank: %d, my number:%d \n",rank, A);

    MPI_Reduce(&A, &somma, 1,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&A, &prodotto, 1,MPI_INT, MPI_PROD, 0, MPI_COMM_WORLD);

    if(rank == 0){
    	printf("Somma= %d, Prodotto=%d", somma, prodotto);
	    printf("tempo_impiegato=%f", MPI_Wtime()-w_time);

    }


MPI_Finalize();


}




