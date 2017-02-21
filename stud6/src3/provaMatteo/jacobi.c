#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include "nrutil.h"



int size, rank, rc;
MPI_Status status; /*Variabile di stato*/
MPI_Request request1, request2, request3, request4, request5, request6, request7, request8, request9, request10;


int main(argc, argv)
int argc;
char *argv[];{
int n = 7;
float s,tau;

float dati[10];

  float a_array[] = {
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    };

    rc=MPI_Init(&argc, &argv);

    if(rc!=MPI_SUCCESS) {
            printf("Inizialization error.\n");
            MPI_ABORT(MPI_COMM_WORLD, &size);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double w_time = MPI_Wtime();
    
if(rank==0) {
    int i, j, nrot;

    float *d = vector(1, n);
    float *r = vector(1, n);

    float **v = matrix(1, n, 1, n);
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
                d[iq] += h; 

                a[ip][iq]=0.0;

                dati[1] =ip ;
                dati[2] =iq;
                dati[3] = 0;
                dati[4]=tau;
                dati[5]=s;

                //PRIMO FOR    
                if( ip-1 > size -1 ){
                    dati[8] = 0;
                    for (j=1;j<size;j++) {
                        dati[0] =j;
                        dati[8] = a[(int)dati[0]][ip];
                        dati[9] = a[(int)dati[0]][iq];
                        MPI_Isend(&dati,10,MPI_FLOAT,j,1, MPI_COMM_WORLD, &request2);
                                       
                    }
                    dati[8] = 1;
                    int contatore=(size);
                    while (contatore<=ip-1) { 
                        MPI_Recv(&dati,10, MPI_FLOAT,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);
                        int invia = status.MPI_SOURCE;

                        a[(int)dati[0]][(int)dati[1]]=dati[6]; //aggiornamento matrice A
                        a[(int)dati[0]][(int)dati[2]]=dati[7];
                        dati[0] = contatore;
                        dati[8] = a[(int)dati[0]][ip];
                        dati[9] = a[(int)dati[0]][iq];
                        MPI_Send(&dati,10,MPI_FLOAT,invia,1, MPI_COMM_WORLD);

                        contatore++;
                    }

                    contatore=1;
                    while(contatore<size){
                        MPI_Recv(&dati,10, MPI_FLOAT,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);
                        a[(int)dati[0]][(int)dati[1]]=dati[6]; //aggiornamento matrice A
                        a[(int)dati[0]][(int)dati[2]]=dati[7];
                        contatore++;
                    }    
                }  
                else{
                    for (j=1;j<=ip-1;j++) {   
                        dati[0] = j; 
                        dati[8] = a[(int)dati[0]][ip];
                        dati[9] = a[(int)dati[0]][iq];
                        MPI_Isend(&dati,10,MPI_FLOAT,j,1, MPI_COMM_WORLD, &request2);
                        MPI_Recv(&dati,10, MPI_FLOAT,j,2,MPI_COMM_WORLD,&status);
                        a[(int)dati[0]][(int)dati[1]]=dati[6]; //aggiornamento matrice A
                        a[(int)dati[0]][(int)dati[2]]=dati[7];
                        //ROTATE(a,j,ip,j,iq)  
                    }
                    
                }
                //FINE PRIMO FOR   
                
                //SECONDO FOR
                int indice=1;
                if( ((iq-1)-(ip+1) +1) > (size-1) )  { //j=ip+1;j<=iq-1;j++; NUM_ITERAZIONI
                    indice=ip+1;
                    for (j=1;j<size;j++) {
                        dati[0] = indice;
                        indice++;
                        dati[8] = a[ip][(int)dati[0]];
                        dati[9] = a[(int)dati[0]][iq];
                        MPI_Isend(&dati,10,MPI_FLOAT,j,1, MPI_COMM_WORLD, &request4);  
                       
                    }

                    int contatore=indice;
                    while (contatore<=iq-1) { 
                        MPI_Recv(&dati,10, MPI_FLOAT,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);
                        int invia = status.MPI_SOURCE;

                        a[(int)dati[1]][(int)dati[0]]=dati[6];
                        a[(int)dati[0]][(int)dati[2]]=dati[7];  
                        dati[0] = contatore;
                        dati[8] = a[ip][(int)dati[0]];
                        dati[9] = a[(int)dati[0]][iq];
                        MPI_Send(&dati,10,MPI_FLOAT,invia,1, MPI_COMM_WORLD);

                        contatore++;
                    }

                    contatore=1;
                    while(contatore<size){
                        MPI_Recv(&dati,10, MPI_FLOAT,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);

                        a[(int)dati[1]][(int)dati[0]]=dati[6];
                        a[(int)dati[0]][(int)dati[2]]=dati[7];  

                        contatore++;
                    }

                }
                else{
                    int contatore=1;
                
                    for (j=ip+1;j<=iq-1;j++) { 
                        
                        dati[0] =j; 
                        dati[8] = a[ip][(int)dati[0]];
                        dati[9] = a[(int)dati[0]][iq];
                        MPI_Isend(&dati,10,MPI_FLOAT,contatore,1, MPI_COMM_WORLD, &request4);
                        MPI_Recv(&dati,10, MPI_FLOAT,contatore,2,MPI_COMM_WORLD,&status);
                        a[(int)dati[1]][(int)dati[0]]=dati[6];
                        a[(int)dati[0]][(int)dati[2]]=dati[7];  
                        contatore++;
                        //ROTATE(a,ip,j,j,iq)
                    }

                  
                }
                //FINE SECONDO FOR
           
                //TERZO FOR
                indice=1;
                if( ((n)-(iq+1) +1) > (size-1) )  { //j=iq+1;j<=n;j++

                    indice=iq+1;
                    for (j=1;j<size;j++) {
                        dati[0] =indice;
                        indice++;
                        dati[8] = a[ip][(int)dati[0]];
                        dati[9] = a[iq][(int)dati[0]];
                        MPI_Send(&dati,10,MPI_FLOAT,j,1, MPI_COMM_WORLD);                         
                    }

                    int contatore=indice; 
                    while (contatore<=n) { 
                        MPI_Recv(&dati,10, MPI_FLOAT,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);
                        int invia = status.MPI_SOURCE;

                        a[(int)dati[1]][(int)dati[0]]=dati[6];
                        a[(int)dati[2]][(int)dati[0]]=dati[7];  
                        dati[0] = contatore;
                        dati[8] = a[ip][(int)dati[0]];
                        dati[9] = a[iq][(int)dati[0]];

                        MPI_Send(&dati,10,MPI_FLOAT,invia,1, MPI_COMM_WORLD);

                        contatore++;
                    }

                    contatore=1;
                    while(contatore<size){
                        MPI_Recv(&dati,10, MPI_FLOAT,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);

                        a[(int)dati[1]][(int)dati[0]]=dati[6];
                        a[(int)dati[2]][(int)dati[0]]=dati[7]; 
                        contatore++;
                    }

                }
                else{
                    int contatore=1;
                    for (j=iq+1;j<=n;j++) { 
                        
                        dati[0] =j; 
                        dati[8] = a[ip][(int)dati[0]];
                        dati[9] = a[iq][(int)dati[0]];
                        MPI_Send(&dati,10,MPI_FLOAT,contatore,1, MPI_COMM_WORLD);
                        MPI_Recv(&dati,10, MPI_FLOAT,contatore,2,MPI_COMM_WORLD,&status);
                        a[(int)dati[1]][(int)dati[0]]=dati[6];
                        a[(int)dati[2]][(int)dati[0]]=dati[7]; 
                        contatore++;   
                        //ROTATE(a,ip,j,iq,j) 
                    }  
                }
                //FINE TERZO FOR


                //QUARTO FOR
                if(n > size-1){
 
                    for (j=1;j<size;j++) {
                        dati[0] =j;
                        dati[8] = v[(int)dati[0]][ip];
                        dati[9] = v[(int)dati[0]][iq];
                        MPI_Send(&dati,10,MPI_FLOAT,j,1, MPI_COMM_WORLD);  
                    }

                    int contatore=(size);
                    while (contatore<=n) { 
                        MPI_Recv(&dati,10, MPI_FLOAT,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);
                        int invia = status.MPI_SOURCE;

                        v[(int)dati[0]][(int)dati[1]]=dati[6];
                        v[(int)dati[0]][(int)dati[2]]=dati[7];
                        dati[0] = contatore;
                        dati[8] = v[(int)dati[0]][ip];
                        dati[9] = v[(int)dati[0]][iq];
                        MPI_Send(&dati,10,MPI_FLOAT,invia,1, MPI_COMM_WORLD);
                        contatore++;  
                    }
                    
                    contatore=1;
                    while(contatore<size){
                        MPI_Recv(&dati,10, MPI_FLOAT,MPI_ANY_SOURCE,2,MPI_COMM_WORLD,&status);

                        v[(int)dati[0]][(int)dati[1]]=dati[6];
                        v[(int)dati[0]][(int)dati[2]]=dati[7];

                        contatore++;
                    }       
                }
                else{
                    for (j=1;j<=n;j++) {    
                        dati[0] = j;
                        dati[8] = v[(int)dati[0]][ip];
                        dati[9] = v[(int)dati[0]][iq];
                        MPI_Send(&dati,10,MPI_FLOAT,j,1, MPI_COMM_WORLD);  
                        MPI_Recv(&dati,10, MPI_FLOAT,j,2,MPI_COMM_WORLD,&status);

                        v[(int)dati[0]][(int)dati[1]]=dati[6];
                        v[(int)dati[0]][(int)dati[2]]=dati[7];
                        
                        //ROTATE(v,j,ip,j,iq)
                    } 
                }
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
    free_matrix(v, 1, n, 1, n);
    free_vector(r, 1, n);
    free_vector(d, 1, n);

    printf("\n");
    printf("tempo impiegato: %f ", tempo_fine);


    //ine
    dati[3]= -1;   
    for(j=1; j<size; j++)
        MPI_Send(&dati,10,MPI_FLOAT,j,1, MPI_COMM_WORLD);
    
    
 } 
    //SLAVE
    else {
    
    float g,h;
    float ricevuti[10];

    float **matrA = matrix(1, n, 1, n);
    float **matrV = matrix(1, n, 1, n);

    
    while(1){
        MPI_Recv(&ricevuti,10, MPI_FLOAT,0,1,MPI_COMM_WORLD,&status);    
        
        if(ricevuti[3] == -1) break;

        if((int)ricevuti[3]==0){ 

            //ROTATE
            g=ricevuti[8]; //ricevuti[0]=j, ricevuti[1]=ip, ricevuti[2]=iq
            h=ricevuti[9]; 
            ricevuti[6]=g-ricevuti[5]*(h+g*ricevuti[4]);
            ricevuti[7]=h+ricevuti[5]*(g-h*ricevuti[4]);

        }
        
        MPI_Send(&ricevuti,10,MPI_FLOAT,0,2, MPI_COMM_WORLD);
 }

}

MPI_Finalize();

return 0;
}
