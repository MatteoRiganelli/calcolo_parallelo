#include <stdio.h>
#include <mpi.h>

int main (int argc,char *argv[]) {
//int argc;
//char *argv[]; {
int size,rank,rc=MPI_INIT (&argc,&argv);
MPI_Status status;
MPI_Request request1, request2, request3, request4,request5;
int err;
MPI_Comm_size(MPI_COMM_WORLD, &size);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

double wtime = MPI_Wtime();

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
        int flag=0;	
       int matrice3[5][5];
	//mando la seconda matrice a tutti
	for(i=1;i<size;i++) {
		MPI_Isend(&matrice2[0][0],20,MPI_INTEGER,i,2,MPI_COMM_WORLD,&request1);
	}
	
	//mando le righe della matrice a 
	for (i=1;i<size;i++) {
		for (j=0;j<4;j++)
			vtemp[j]=matrice1[i-1][j];
		vtemp[4]=i-1;
		
		MPI_Isend(&vtemp,5,MPI_INTEGER,i,1,MPI_COMM_WORLD,&request5); //modulo su i e size. Send NON bloccante.
                MPI_Test(&request5,&flag,&status);
                while (flag==0) {
                  MPI_Test(&request5,&flag,&status); 
                }
                flag=0;
	}
	
	int contatore=(size);
	int rcontatore=0;
	while (contatore<=5) {
		
		MPI_Recv(&vtemp,5,MPI_INTEGER,MPI_ANY_SOURCE,3,MPI_COMM_WORLD, &status);
	
		for (i=0;i<4;i++)
			matrice3[vtemp[4]][i]=vtemp[i];	
			for(i=0;i<4;i++)				
				vtemp[i]=matrice1[contatore-1][i];
			vtemp[4]=contatore-1;
			MPI_Isend(&vtemp,5,MPI_INTEGER,status.MPI_SOURCE,1,MPI_COMM_WORLD,&request5); 
			contatore++;
		} 
		contatore=1;
		
		while(contatore<size){
			MPI_Recv(&vtemp,5,MPI_INTEGER,MPI_ANY_SOURCE,3,MPI_COMM_WORLD, &status);
			
			for (i=0;i<4;i++)
				matrice3[vtemp[4]][i]=vtemp[i];
			vtemp[4]=-1;
			MPI_Isend(&vtemp,5,MPI_INTEGER,status.MPI_SOURCE,1,MPI_COMM_WORLD,&request4);
			contatore++;
		}
		

	for(i=0;i<5;i++)
	 	printf("%d,%d,%d,%d \n",matrice3[i][0],matrice3[i][1],matrice3[i][2],matrice3[i][3]);
printf("Tempo: %f \n", MPI_Wtime()-wtime);


//SLAVE
}else{
	int i,j;
	int matrice [4][5]; 

	MPI_Irecv(&matrice[0][0],20,MPI_INTEGER,0,2,MPI_COMM_WORLD,&request2);
	
	//MPI_Wait(&request1, &status,&err);

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
				MPI_Isend(&vet[0],5,MPI_INTEGER,0,3,MPI_COMM_WORLD, &request4);
				
		} 

}

MPI_Finalize ();

}
