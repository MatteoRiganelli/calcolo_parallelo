#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include "nrutil.h"


int size, rank, rc;
MPI_Status status; /*Variabile di stato*/


int main(argc, argv)
int argc;
char *argv[];{
    int n = 2;
    float s,tau;

    double w_time = MPI_Wtime();
    

    rc=MPI_Init(&argc, &argv);

    if(rc!=MPI_SUCCESS) {
            printf("Inizialization error.\n");
            MPI_ABORT(MPI_COMM_WORLD, &size);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    

if(rank==0) {
     float **a = matrix(1, n, 1, n);
     int ip, iq;
     for(ip=1;ip<=n;ip++) {  
        for (iq=1;iq<=n;iq++) a[ip][iq]=0.0; 
        a[ip][ip]=1.0;
    }
    int j;
    for(j=1; j<size; j++)
        MPI_Send(&a[1][1],4,MPI_FLOAT,j,0, MPI_COMM_WORLD);  
        //MPI_Send(&v[0][0],16,MPI_FLOAT,j,0, MPI_COMM_WORLD);  
        //MPI_Send(&t,2,MPI_FLOAT,j,0, MPI_COMM_WORLD);
}


 

   else {
     float **matrice = matrix(1, n, 1, n);
    MPI_Recv(&matrice[1][1],4, MPI_FLOAT,0,0,MPI_COMM_WORLD,&status);
    printf("\n");
    printf("Matrix:\n");
    int i, j;
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            printf("  %f", matrice[i][j]);
        }
        printf("\n");
    }
    printf("\n");
 }

MPI_Finalize();

return 0;

}

