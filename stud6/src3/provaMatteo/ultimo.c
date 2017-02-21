#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include "nrutil.h"
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
a[k][l]=h+s*(g-h*tau);


int size, rank, rc;
MPI_Status status; /*Variabile di stato*/
MPI_Request request1, request2, request3, request4, request5, request6, request7, request8, request9, request10;


int main(argc, argv)
int argc;
char *argv[];{
int n =1700;
float s,tau;
int i, j, nrot;

float dati[10];

    rc=MPI_Init(&argc, &argv);

    if(rc!=MPI_SUCCESS) {
            printf("Inizialization error.\n");
            MPI_ABORT(MPI_COMM_WORLD, &size);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double w_time = MPI_Wtime();
  

if(rank==0) {

    //printf("SIZE =%d ", size );

    /*
       We need a Numerical Recipes-style "matrix" to use
       as an argument to call the jacobi function.
       Here are the matrix values stored in row-major
       order in a 1-D array. This array name will be used
       as an argument for the convert_matrix function.
     */

    float *d = vector(1, n);
    float *r = vector(1, n);

    //float **v = matrix(1, n, 1, n);
    //float **a = convert_matrix(&a_array[0], 1, n, 1, n);
    float **a = matrix(1, n, 1, n);
    
    int axx=1;
    int ayy=1;
    for(axx=1; axx<=n; axx++)
        for(ayy=1; ayy<=n; ayy++)
            a[axx][ayy]=1.0;
    
 /*Jacobi*/
   
int ip,iq;
float tresh,theta,t,sm,h,g,c,*b,*z;

b=vector(1,n); z=vector(1,n);


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
                d[iq] += h; 

                a[ip][iq]=0.0;

                dati[1] =ip ;
                dati[2] =iq;
                dati[3] = 0;
                dati[4]=tau;
                dati[5]=s;

                
//PRIMO FOR    
                
                
                
                    for (j=1;j<=ip-1;j++) {   
                        ROTATE(a,j,ip,j,iq)  
                    }
                    
                
                
                //FINE PRIMO FOR   
                
                
                //SECONDO FOR
                
                    for (j=ip+1;j<=iq-1;j++) { 
                        

                        ROTATE(a,ip,j,j,iq)

                    }

                  
                
                //FINE SECONDO FOR
           
                //TERZO FOR
                    for (j=iq+1;j<=n;j++) { 
                        ROTATE(a,ip,j,iq,j) 
                    }
                    
                
                
                //FINE TERZO FOR


                //QUARTO FOR
                    for(j=1; j< size; j++){
                        MPI_Send(&dati,10,MPI_FLOAT,j,1, MPI_COMM_WORLD);  
                        //MPI_Recv(&v[(int)indiceVec[0]][0],2*n, MPI_FLOAT,j,3,MPI_COMM_WORLD,&status);
                    }        
                //VEDI SLAVE
                
                //FINE QUARTO FOR            
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
    double tempo_fine = MPI_Wtime() - w_time;
    
    printf("Eigenvalues:\n");
    for (i = 1; i <= n; i++) {
        printf("  Number %d:%+14.6e\n", i, d[i]);
    }
    printf("\n");

    
    

    printf("Number of Jacobi rotations = %d\n", nrot);
    free_matrix(a, 1, n, 1, n);
    free_vector(r, 1, n);
    free_vector(d, 1, n);

    printf("\n");
    printf("tempo impiegato: %f ", tempo_fine);


    dati[3]= -1;   
    for(j=1; j<size; j++)
        MPI_Send(&dati,10,MPI_FLOAT,j,1, MPI_COMM_WORLD);
    
    
} 
   else {

   int chunk = n / (size-1);
   int dispari = n % (size-1);
 
    float g,h;
    float ricevuti[10];
    int indici[(size-1)*2];
    float **v2 = matrix(1, n, 1, n);
    double tempo_fine = 0;

    for(i=1;i<=n;i++) {  
        for (j=1;j<=n;j++) v2[i][j]=0.0; 
        v2[i][i]=1.0;
    }


    int indiceIni = 1;
    int indiceFin = chunk;
   
    int indice = 0;

    for(j=1; j<size; j++){
       
        indici[indice] = indiceIni;
        indici[indice+1] = indiceFin;

        indiceIni = indici[indice+1] + 1;

        if(dispari == 0){
            indiceFin = indici[indice+1] + chunk;
        }
        else{

            if(j == size -1) indiceFin = n;
            else indiceFin = indici[indice+1] + chunk;
        }
        
        indice = indice +2 ; 
    }

    while(1){
        
        MPI_Recv(&ricevuti,10, MPI_FLOAT,0,1,MPI_COMM_WORLD,&status); 

        if(ricevuti[3] == -1){ 
                  
                //printf("\nRANGO=%d\n", rank );
                printf("Eigenvectors:\n");
                for (i = 1; i <= n; i++) {
                    printf("  Number %d:", i);
                    for (j = 1; j <= n; j++) {
                        if(v2[j][i] != 0 && v2[j][i]!=1)
                            printf("%+14.6e", v2[j][i]);
                    }
                    printf("\n");
                }
                printf("\n");
                free_matrix(v2, 1, n, 1, n);
            break;
        }
        
            int indice1= (2*rank)-2;
            int indice2=(2*rank)-1;
            //printf("RANGO=%d - indice1=%d, indice2=%d\n", rank, indice1, indice2);
            for(j=indici[indice1]; j<= indici[indice2]; j++){
                    //ROTATE(v2,j, (int)dati[1], j, (int)dati[2]);
                g=v2[j][(int)ricevuti[1]]; //ricevuti[0]=j, ricevuti[1]=ip, ricevuti[2]=iq
                h=v2[j][(int)ricevuti[2]]; 
                ricevuti[6]=g-ricevuti[5]*(h+g*ricevuti[4]);
                ricevuti[7]=h+ricevuti[5]*(g-h*ricevuti[4]);
                v2[j][(int)ricevuti[1]] =  ricevuti[6];
                v2[j][(int)ricevuti[2]] =  ricevuti[7];
                

            }
        
            //MPI_Send(&v2[(int)indici[2*rank-1]][0], 2*n, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);    
        
 }

 tempo_fine = MPI_Wtime() - w_time;
printf("tempo slave - %d = %f\n", rank, tempo_fine);



}

MPI_Finalize();

return 0;
}
