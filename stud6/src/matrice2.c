#include <mpi.h>
#include <stdio.h>


int main(argc, argv)
int argc;
char *argv[];{
	int size,rank, rc, i, j, k, t, pr;
	MPI_Status status;
	int matA[4][4] = {2,2,2,4,
		         5,9,8,5,
			 4,7,5,8,
			 9,0,2,2};
int matB[4][4] = {1,0,0,0,
		  0,1,0,0,
		  0,0,1,0,
		  0,0,0,1};
/*
int matB[4][4] = {4,7,8,4,
		  5,6,3,2,
		  1,5,9,6,
		  2,0,4,7};
*/					  
	int matC[4][4];
	int matR[4][4];
	int array1[4];
	int array2[4];
	
	rc = MPI_INIT(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	if(rank == 0){
		for(j=1; j<5; j++){
			MPI_Send(&matB[0][0], 16, MPI_INTEGER, j, 0, MPI_COMM_WORLD); 
			MPI_Send(&matA[j-1][0], 4, MPI_INTEGER, j, 0, MPI_COMM_WORLD);

		}
		
		//receive del master (rank 0)
		for(j=1; j<5; j++){	
			MPI_Recv(&matR[j-1][0], 4, MPI_INTEGER, j, 0, MPI_COMM_WORLD, &status);
		  
		}
		printf("RISULTATO \n");
		for(i=0; i<4; i++){
		printf("\n");
			for(j=0; j<4; j++){
				printf("%d ", matR[i][j]);
			}
		}
		
		
	}else{
		for(t=1; t<5; t++){
			if(t == rank){
				MPI_Recv(&matC[0][0], 16, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&array1[0],4 , MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &status);
				/*
				printf("ARRAY1 proc %d -> ", t); 
				for(pr=0; pr < 4; pr++)
					printf(" %d", array1[pr]);
				printf("\n");
				*/
				
				for(i=0; i<4; i++){
					array2[i] = 0;
					for(j=0; j<4; j++){
					array2[i]+=array1[j]*matC[j][i];
					}
					
				}
				/*
				printf("ARRAY2 proc %d -> ", t); 
					for(pr=0; pr < 4; pr++)
						printf(" %d", array2[pr]);
				printf("\n");
				*/
				MPI_Send(&array2[0], 4, MPI_INTEGER,0,0,MPI_COMM_WORLD);  
			}
		}
	
	}
	
	MPI_Finalize();
	
}
	
