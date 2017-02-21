#include <mpi.h>
#include <stdio.h>

int main(argc, argv)

	int argc;
	char *argv[];   {

	int size, rank, rc, i,j,dest, x;       
            
        MPI_Status status; /*Variabile di stato*/

         int n=4,m=4; /*Righe e colonne matrice*/
         int matrC[n][m]; /*La Matrice che ricevono i task*/
         int ris[n][m];  /*La matrice risultato*/
         int riga[n];  /*la riga ricevuta dal task*/
         int rigaTask[n];   /*La riga della matrice ris che viene inviata al master*/
         /* Inizializza riga task vuota*/
         for (i=0;i<n;i++) {
           rigaTask[i]=0;
          }

	rc=MPI_Init(&argc, &argv);

	if(rc!=MPI_SUCCESS) {
		printf("Inizialization error.\n");
		MPI_ABORT(MPI_COMM_WORLD, &size);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
        /*Codice che eseguirà Il master*/

          /*Inizializza e visualizza Matrice 1*/

          int matr1[4][4];

          /*Inizializza e visualizza Matrice 2*/

         int matr2[4][4];


if(rank==0) {

  for(i=0;i<4;i++) {
   for(j=0;j<4;j++) {
     matr1[i][j]=1;
     matr2[i][j]=1;
   }
  }

}

  /*Invia la matrice B ai task figli e poi una riga ciascuno*/

   MPI_Bcast(&matr2,16,MPI_INT,0,MPI_COMM_WORLD);


   MPI_Scatter(&matr1,4,MPI_INT,&riga,4,MPI_INT,0,MPI_COMM_WORLD); 

      
       for(i=0;i<n;i++) {
         for (x=0; x<m;x++){
           rigaTask[i]+= riga[x]*matr2[x][i];
           printf(" %d ",rigaTask[i]);  
         }      
        printf("\n");
       }
  

printf("La mia riga %d \n",rank);
/*
MPI_Gather(&ris[rank],4,MPI_INT,&ris,4,MPI_INT,0,MPI_COMM_WORLD);


if (rank==0) {    
       
      for (i=0;i<n;i++) {
        for (j=0;j<m;j++) {
          printf(" %d ",ris[i][j]);
        }
        printf("\n");

      }
}



printf("Sono il task %d \n  ",rank);
for (i=0;i<n;i++) {
        for (j=0;j<m;j++) {
          printf(" %d ",matr2[i][j]);
        }
	printf("\n");

      }

printf("\n ");
*/
MPI_Finalize();

}
