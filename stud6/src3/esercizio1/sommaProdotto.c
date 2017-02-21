#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

int main (int argc,char *argv[]) {


int size,rank,rc=MPI_INIT (&argc,&argv);
double tempo=MPI_Wtime();
MPI_Status status;
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Request request;
	
	srand(time(NULL) + rank);
	int sommas=(rand()%100) +1;	
	long int prodottos=(rand()%10) +1;
	long int prodottor;
	int sommar;
	int passo=4;
	//printf("valore: %d\n",a);
   	
	if (rank%2==0){
		MPI_Recv(&sommar,1,MPI_INTEGER,rank+1,1,MPI_COMM_WORLD, &status);
		MPI_Recv(&prodottor,1,MPI_INTEGER,rank+1,1,MPI_COMM_WORLD, &status);
		sommas=sommas+sommar;
		prodottos=prodottos*prodottor;
		printf("ricevente: %d somma: %d prodotto: %d\n", rank,sommas,prodottos);
}else{
		MPI_Send(&sommas,1,MPI_INTEGER,rank-1,1,MPI_COMM_WORLD);
		MPI_Send(&prodottos,1,MPI_INTEGER,rank-1,1,MPI_COMM_WORLD);
}
	while(passo<=size){
		if (rank%passo==0){//ricevente
			MPI_Recv(&sommar,1,MPI_INTEGER,rank+passo/2,1,MPI_COMM_WORLD, &status);
			MPI_Recv(&prodottor,1,MPI_INTEGER,rank+passo/2,1,MPI_COMM_WORLD, &status);
			sommas=sommas+sommar;
			prodottos=prodottos*prodottor;
			printf("ricevente: %d somma: %d prodotto: %d\n", rank,sommas,prodottos);
		}else if (rank%2==0 && rank%passo!=0){
			MPI_Send(&sommas,1,MPI_INTEGER,rank-passo/2,1,MPI_COMM_WORLD); 
			MPI_Send(&prodottos,1,MPI_INTEGER,rank-passo/2,1,MPI_COMM_WORLD);
		}
		passo=passo*2;
	}
	
	if(rank==0){
		double tempo2=MPI_Wtime();
		printf("la somma è: %d prodotto è: %d  e il tempo è:%f\n",sommas,prodottos,(tempo2-tempo));
	}
	MPI_Finalize();
}
