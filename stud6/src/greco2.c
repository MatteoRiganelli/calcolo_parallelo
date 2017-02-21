#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int main(argc, argv)

        int argc;
        char *argv[];   {

        double w_time = MPI_Wtime();

        int size, rank, rc,dest;
        int N=10000;
        int punti_cerchio=0;

        MPI_Status status; /*Variabile di stato*/

        rc=MPI_Init(&argc, &argv);

        if(rc!=MPI_SUCCESS) {
                printf("Inizialization error.\n");
                MPI_ABORT(MPI_COMM_WORLD, &size);
        }

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

       int num=N/7;

       int j;
       double cord1,cord2;

    for(j=1;j<num;j++) {
      cord1=rand()/(double)RAND_MAX;
      cord2=rand()/(double)RAND_MAX;

      if(cord1*cord1+cord2*cord2 <= 1 )
        ++punti_cerchio;
      }

if (rank==0) {
     int j;
     int ris[6];

   for(j=1;j<7;j++) {
      MPI_Recv(&ris[j-1],1, MPI_INTEGER,j,0, MPI_COMM_WORLD, &status);
   }

   int somma=0;

   for(j=0;j<6;j++) {
    // printf(" %d ",ris[j]);
     somma+=ris[j];
  }

double finale= ((double)(somma+punti_cerchio)/N)*4;

printf(" %f\n ", finale );

printf("tempo impiegato: %f ", MPI_Wtime() - w_time);

}

 else  {
//Worker
    
   // printf(" %f \n",cord1);
   // printf("my rank  %d \n",rank);

 MPI_Send(&punti_cerchio,1,MPI_INTEGER,0,0,MPI_COMM_WORLD);


}

MPI_Finalize();


}


