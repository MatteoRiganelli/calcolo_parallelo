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

if (rank==0) {

          /*Inizializza e visualizza Matrice 1*/

          int matr1[4][4] = {1,1,2,4,
                            5,9,8,5,
                            4,7,5,8,
                            9,0,2,2};

          /*Inizializza e visualizza Matrice 2*/

         int matr2[4][4] = {1,0,0,0,
                           0,1,0,0,
                           0,0,1,0,
                           0,0,0,1};
     
  /*Invia la matrice B ai task figli e poi una riga ciascuno*/

       for(dest=1;dest < size; dest++) {
          MPI_Send(&matr2[0][0],16,MPI_INTEGER,dest,0, MPI_COMM_WORLD);         
          MPI_Send(&matr1[dest-1][0],4,MPI_INTEGER,dest,0,MPI_COMM_WORLD);
         
       } 
     
      /*Ricevi le righe */
       for(dest=1;dest<size;dest++ ) {
         MPI_Recv(&ris[dest-1][0],4,MPI_INTEGER,dest,0,MPI_COMM_WORLD,&status);
       } 

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
      
       /*Ricevi la riga e ricevi la matrice B*/
       MPI_Recv(&matrC[0][0], 16, MPI_INTEGER,0, 0, MPI_COMM_WORLD, &status);
       MPI_Recv(&riga[0],4, MPI_INTEGER,0,0,MPI_COMM_WORLD,&status);
          
       /*Calcola i prodotti riga per colonna e restituisci al Master*/
       for(i=0;i<n;i++) {
         for (x=0; x<m;x++){
           rigaTask[i]+= riga[x]*matrC[x][i];  
         }      
       }
       
     MPI_Send(&rigaTask[0],4,MPI_INT,0,0,MPI_COMM_WORLD);
       
}                

MPI_Finalize();

}
