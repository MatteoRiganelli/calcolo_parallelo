#include <mpi.h>
#include <stdio.h>

int main(argc, argv)

        int argc;
        char *argv[];   {

        int size, rank, rc, i,j,dest, x;

        MPI_Status status; /*Variabile di stato*/

         int n=8,m=8; /*Righe e colonne matrice*/
         int matrC[n][m]; /*La Matrice che ricevono i task*/
         int ris[n][m];  /*La matrice risultato*/
         int riga[n];  /*la riga ricevuta dal task*/
         int rigaTask[n];   /*La riga della matrice ris che viene inviata al master*/
         /* Inizializza riga task vuota*/
         for (i=0;i<n;i++) {
           rigaTask[i]=0;
          }

        int tag1=1,tag2=2, tag3;

        int invia;

        rc=MPI_Init(&argc, &argv);

        if(rc!=MPI_SUCCESS) {
               	printf("Inizialization error.\n");
                MPI_ABORT(MPI_COMM_WORLD, &size);
        }

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

	/*Codice che eseguirà Il master*/
if (rank==0) {
           int tag3=1;
          /*Inizializza e visualizza Matrice 1*/

          int matr1[8][8] = {1,1,2,4,1,1,2,4,
                            5,9,8,5,5,9,8,5,
                            4,7,5,8,4,7,5,8,
                            9,0,2,2,9,0,2,2,
                            1,1,2,4,1,1,2,4,
                            5,9,8,5,5,9,8,5,
                            4,7,5,8,4,7,5,8,
                            9,0,2,2,9,0,2,2};

          /*Inizializza e visualizza Matrice 2*/

         int matr2[8][8] = {1,0,0,0,1,0,0,0,
                           0,1,0,0,0,1,0,0,
                           0,0,1,0,0,0,1,0,
                           0,0,0,1,0,0,0,1,
                           1,0,0,0,1,0,0,0,
                           0,1,0,0,0,1,0,0,
                           0,0,1,0,0,0,1,0,
                           0,0,0,1,0,0,0,1,};

/*Invia la matrice B ai task figli una volta sola e poi una riga ciascuno delle prime 4*/

       for(dest=1;dest < size; dest++) {
          MPI_Send(&matr2[0][0],64,MPI_INTEGER,dest,0, MPI_COMM_WORLD);
       }

       for(dest=1;dest < size; dest++) {
          MPI_Send(&matr1[dest-1][0],8,MPI_INTEGER,dest,i,MPI_COMM_WORLD);
       }

/*Invia le righe rimanenti e aspetta la risposta*/
      for(dest=size;dest<m;dest++ ) {
         MPI_Recv(&ris[dest-1][0],8,MPI_INTEGER,MPI_ANY_SOURCE,tag3,MPI_COMM_WORLD,&status);
         invia= status.MPI_SOURCE;
         MPI_Send(&matr1[dest-1][0],8,MPI_INTEGER, invia,i,MPI_COMM_WORLD);

       }

         MPI_Send(&matr1[0][0],8,MPI_INTEGER, dest,-1,MPI_COMM_WORLD);

       /*mostra*/
      for (i=0;i<n;i++) {
        for (j=0;j<m;j++) {
          printf(" %d ",ris[i][j]);
        }
	printf("\n");
    }
}
 else  {

   /*I processi figli*/
    
   /*Ricevi la matrice B una volta sola */
   MPI_Recv(&matrC[0][0], 64, MPI_INTEGER,0, tag1, MPI_COMM_WORLD, &status);
   
   while(1) {
       
       /*Ricevi la riga e ricevi la matrice B*/
       MPI_Recv(&riga[0],8, MPI_INTEGER,0,tag3,MPI_COMM_WORLD,&status);
       
      // if (tag3==0) {
       /*Calcola i prodotti riga per colonna e restituisci al Master*/
       for(i=0;i<n;i++) {
         for (x=0; x<m;x++){
           rigaTask[i]+= riga[x]*matrC[x][i];
         }
       }
          MPI_Send(&rigaTask[0],8,MPI_INT,0,0,MPI_COMM_WORLD);
      // } else if (tag3==-1) {
        //  break;
      // }
   }
}

MPI_Finalize();

}
