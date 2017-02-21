#include <stdio.h>
#include <mpi.h>

int main (int argc,char *argv[]) {
//int argc;
//char *argv[]; {
int size,rank,rc=MPI_INIT (&argc,&argv);
MPI_Status status;
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

if (rank==0) {
	int i,j;
	int vtemp[5];

	int matrice1 [5][4]={{1,1,1,1},
		    			 {2,2,2,2},
		     			 {4,4,4,4},
		     			 {3,3,3,3},
		     			 {4,4,4,4}};

	int matrice2 [4][5]= {{1,1,1,1,1},
	             		  {1,1,1,1,1},
		     			  {1,1,1,1,1},
		     			  {1,1,1,1,1}};
	int matrice3[5][5];
	//mando la seconda matrice a tutti
	for(i=1;i<size;i++) {
		MPI_Send(&matrice2[0][0],20,MPI_INTEGER,i,2,MPI_COMM_WORLD);
	}
	
	//mando le righe della matrice a 
	for (i=1;i<size;i++) {
		for (j=0;j<4;j++)
			vtemp[j]=matrice1[i-1][j];
		vtemp[4]=i-1;
		
		MPI_Send(&vtemp,5,MPI_INTEGER,i,1,MPI_COMM_WORLD); //modulo su i e size. Send NON bloccante.
	}
	
	int contatore=(size);
	printf("contatore= %d\n", contatore);
	int rcontatore=0;
	while (contatore<=5) {
		
		MPI_Recv(&vtemp,5,MPI_INTEGER,MPI_ANY_SOURCE,3,MPI_COMM_WORLD, &status);
	
		for (i=0;i<4;i++)
			matrice3[vtemp[4]][i]=vtemp[i];	
			for(i=0;i<4;i++)				
				vtemp[i]=matrice1[contatore-1][i];
			vtemp[4]=contatore-1;
			MPI_Send(&vtemp,5,MPI_INTEGER,status.MPI_SOURCE,1,MPI_COMM_WORLD); 
			contatore++;
	} 
		contatore=1;
		
		while(contatore<size){
			MPI_Recv(&vtemp,5,MPI_INTEGER,MPI_ANY_SOURCE,3,MPI_COMM_WORLD, &status);
			
			for (i=0;i<4;i++)
				matrice3[vtemp[4]][i]=vtemp[i];
			vtemp[4]=-1;
			MPI_Send(&vtemp,5,MPI_INTEGER,status.MPI_SOURCE,1,MPI_COMM_WORLD);
			contatore++;
		}
		

	for(i=0;i<5;i++)
	 	printf("%d,%d,%d,%d \n",matrice3[i][0],matrice3[i][1],matrice3[i][2],matrice3[i][3]);


//SLAVE
}else{
	int i,j;
	int matrice [4][5]; 
	MPI_Recv(&matrice[0][0],20,MPI_INTEGER,0,2,MPI_COMM_WORLD, &status);
	
		
	while(1){
		int vet[5]={0,0,0,0,0};
		int vet2[5];
		MPI_Recv(&vet2[0],5,MPI_INTEGER,0,1,MPI_COMM_WORLD, &status);
		vet[4]=vet2[4];
		if(vet2[4]==-1)
			break;
		for(i=0;i<4;i++)
			for(j=0;j<4;j++) 
				vet[i]= vet[i]+(vet2[j]*matrice[j][i]);
							
				MPI_Send(&vet[0],5,MPI_INTEGER,0,3,MPI_COMM_WORLD);
			
		} 
}

MPI_Finalize ();

}
