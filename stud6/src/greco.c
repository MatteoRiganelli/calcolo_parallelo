#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(argc, argv)

        int argc;
        char *argv[];   {

	double w_time = MPI_Wtime();

        int size, rank, rc,dest;
        int N=1000000000;
        int punti_cerchio=0;

        MPI_Status status; /*Variabile di stato*/

        rc=MPI_Init(&argc, &argv);

        if(rc!=MPI_SUCCESS) {
               	printf("Inizialization error.\n");
                MPI_ABORT(MPI_COMM_WORLD, &size);
        }

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

       int num=N/5;
	/*Codice che eseguirà Il master*/

if (rank==0) {
     int j;  
     int ris[5];    

   for(j=1;j<size;j++) {
     MPI_Recv(&ris[j-1],1, MPI_INTEGER,j,0, MPI_COMM_WORLD, &status);  
   }

   int somma=0;

   for(j=0;j<5;j++) {
    // printf(" %d ",ris[j]);
     somma+=ris[j];
  }

double finale= ((double)somma/N)*4; 

printf(" %f\n ", finale );

printf("tempo impiegato: %f ", MPI_Wtime() - w_time);

}

 else  {
//Worker
    int j;
    double cord1,cord2;

    for(j=1;j<num;j++) {
      cord1=rand()/(double)RAND_MAX;
      cord2=rand()/(double)RAND_MAX;
    
      if(cord1*cord1+cord2*cord2 < 1 )
        ++punti_cerchio;
      }

   // printf(" %f \n",cord1);
    printf("my rank  %d \n",rank);

    MPI_Send(&punti_cerchio,1,MPI_INTEGER,0,0,MPI_COMM_WORLD);

   /*I processi figli*/
         //MPI_Recv(&matrC[0][0], 64, MPI_INTEGER,0, tag1, MPI_COMM_WORLD, &status);
   
}

MPI_Finalize();


}




