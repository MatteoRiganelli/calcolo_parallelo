#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include "nrutil.h"
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
a[k][l]=h+s*(g-h*tau);


int size, rank, rc;
MPI_Status status; /*Variabile di stato*/

int main(argc, argv)
int argc;
char *argv[];{
    int n = 2;
    float s,tau;

 float a_array[] = {
         5.0, -2.0, 
        -2.0,  2.0

    };

    double w_time = MPI_Wtime();
    

    rc=MPI_Init(&argc, &argv);

    if(rc!=MPI_SUCCESS) {
            printf("Inizialization error.\n");
            MPI_ABORT(MPI_COMM_WORLD, &size);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    

if(rank==0) {
    int i, j, nrot;

    /*
       We need a Numerical Recipes-style "matrix" to use
       as an argument to call the jacobi function.
       Here are the matrix values stored in row-major
       order in a 1-D array. This array name will be used
       as an argument for the convert_matrix function.
     */

   
    float *d = vector(1, n);
    float *r = vector(1, n);
    float v[n][n];
    float **a = convert_matrix(&a_array[0], 1, n, 1, n);
    /*

    printf("Matrix:\n");
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            printf("  %+13.6e", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
*   /
    //jacobi(a, n, d, v, &nrot);
    
     /*Jacobi*/
   
int ip,iq;
float tresh,theta,t,sm,h,g,c,*b,*z;

b=vector(1,n); z=vector(1,n);
for(ip=1;ip<=n;ip++) {  
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0; 
    v[ip][ip]=1.0;
}

for(ip=1;ip<=n;ip++) { 
    b[ip]=d[ip]=a[ip][ip]; 
    z[ip]=0.0;
} 
nrot=0;
for(i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) {
        for (iq=ip+1;iq<=n;iq++) sm += fabs(a[ip][iq]);
    } 
    if(sm==0.0){
        free_vector(z,1,n); 
        free_vector(b,1,n); 
        break;
    }
    if (i < 4)
        tresh=0.2*sm/(n*n); 
    else
        tresh=0.0;
    for (ip=1;ip<=n-1;ip++) {
        for (iq=ip+1;iq<=n;iq++) {
            g=100.0*fabs(a[ip][iq]);
            if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip]) && (float)(fabs(d[iq])+g) == (float)fabs(d[iq])) 
                a[ip][iq]=0.0;
            else if (fabs(a[ip][iq]) > tresh) { 
                h=d[iq]-d[ip];
                if ((float)(fabs(h)+g) == (float)fabs(h)) 
                    t=(a[ip][iq])/h;
                else{
                    theta=0.5*h/(a[ip][iq]); 
                    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                    if (theta < 0.0) t = -t;
                }
                c=1.0/sqrt(1+t*t); 
                s=t*c;
                tau=s/(1.0+c); 
                h=t*a[ip][iq];
                z[ip] -= h;
                z[iq] += h;
                d[ip] -= h;
                d[iq] += h; a[ip][iq]=0.0;

                                                         

                for (j=1;j<=ip-1;j++) {
                    ROTATE(a,j,ip,j,iq)
                }
                for (j=ip+1;j<=iq-1;j++) { 
                    ROTATE(a,ip,j,j,iq)
                }
                for (j=iq+1;j<=n;j++) {
                    ROTATE(a,ip,j,iq,j) 
                }
                                
                                
                                
                int dati[3]; //array che memorizza gli indici da inviare
                float t[2]; //array che memorizza tau, s
                float mettereV[2]; //dove mettere i risultati nella matrice v
                //inizializzo t
                t[0]=tau; 
                t[1]=s;
                for(j=1; j<size; j++){
                    MPI_Send(&v[0][0],4,MPI_FLOAT,j,0, MPI_COMM_WORLD); //invio v a tutti gli slave
                    MPI_Send(&t,2,MPI_FLOAT,j,0, MPI_COMM_WORLD); //invio tao e s a tutti gli slave
                }
                //per ogni 
                for (j=1;j<=n;j++) {
                                         dati[0] =j ;
                                         dati[1] =ip ;
                                         dati[2] =iq;
                                         MPI_Send(&dati,3,MPI_INTEGER,j,0, MPI_COMM_WORLD);
                                         MPI_Recv(&dati,3, MPI_INTEGER,j,0,MPI_COMM_WORLD,&status);
                                         MPI_Recv(&mettereV,2, MPI_FLOAT,j,0,MPI_COMM_WORLD,&status);
                                         v[dati[0]][dati[1]]=mettereV[0];
                                         v[dati[0]][dati[2]]=mettereV[1];
                                         //ROTATE(v,j,ip,j,iq)
                                   
                }


                //MPI_Recv(&v[0][0],4, MPI_INT,0,0,MPI_COMM_WORLD,&status);
 

                ++nrot; 
            }
        }
    } 
    for (ip=1;ip<=n;ip++){ 
        b[ip] += z[ip]; 
        d[ip]=b[ip]; 
        z[ip]=0.0;
                
    }

}
 
    printf("Eigenvalues:\n");
    for (i = 1; i <= n; i++) {
        printf("  Number %d:%+14.6e\n", i, d[i]);
    }
    printf("\n");

    printf("Eigenvectors:\n");
    for (i = 1; i <= n; i++) {
        printf("  Number %d:", i);
        for (j = 1; j <= n; j++) {
            printf("%+14.6e", v[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    printf("Number of Jacobi rotations = %d\n", nrot);

    free_convert_matrix(a, 1, n, 1, n);
    //free_matrix(v, 1, n, 1, n);
    free_vector(r, 1, n);
    free_vector(d, 1, n);

    printf("\n");
    printf("tempo impiegato: %f ", MPI_Wtime() - w_time);



 } 

   else {

     int i,j;
     float g,h;
   
       int ricevuti[3];
       float matrV[n][n];
       float t[2];
        MPI_Recv(&matrV[0][0],4, MPI_FLOAT,0,0,MPI_COMM_WORLD,&status);
        MPI_Recv(&t,2, MPI_FLOAT,0,0,MPI_COMM_WORLD,&status);

        printf("Matrix:\n");
        for (i = 0; i <n; i++) {
            for (j = 0; j < n; j++) {
                printf("  %f", matrV[i][j]);
            }
            printf("\n");
        }
        printf("\n");
        MPI_Recv(&ricevuti,3, MPI_INT,0,0,MPI_COMM_WORLD,&status);

        //ROTATE(j, ip, j, iq)
        //ricevuti[0]=j;
        //ricevuti[1]=ip;
        //ricevuti[2]=iq;
        //ROTATE(matrV,ricevuti[0],ricevuti[1],ricevuti[0],ricevuti[2])
        g=matrV[ricevuti[0]][ricevuti[1]];
        h=matrV[ricevuti[0]][ricevuti[2]];
        matrV[ricevuti[0]][ricevuti[1]]=g-t[1]*(h+g*t[0]);
        matrV[ricevuti[0]][ricevuti[2]]=h+t[1]*(g-h*t[0]);
        printf("a[i][j]= %+14.6e ", matrV[ricevuti[0]][ricevuti[1]]);
        printf("a[k][l]= %+14.6e\n", matrV[ricevuti[0]][ricevuti[2]]);
        MPI_Send(&ricevuti,3,MPI_INTEGER,0,0, MPI_COMM_WORLD);

        float ris[2];
        ris[0]=matrV[ricevuti[0]][ricevuti[1]];
        ris[1]=matrV[ricevuti[0]][ricevuti[2]];

        printf("ris0= %+14.6e\n, ris1= %+14.6e\n\n", ris[0], ris[1] );

        MPI_Send(&ris,2,MPI_FLOAT,0,0, MPI_COMM_WORLD);
 }

MPI_Finalize();

return 0;
}
