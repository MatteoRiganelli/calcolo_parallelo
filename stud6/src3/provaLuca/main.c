#include <math.h>
#include "nrutil.h"
#include <mpi.h>
#include <stdio>

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
a[k][l]=h+s*(g-h*tau);

int j,iq,ip,i;
float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;


int main(argc, argv)
int argc;
char *argv[];{
   // double w_time = MPI_Wtime();
    
     int i, j, nrot;

    int n = 2;
    /*
       We need a Numerical Recipes-style "matrix" to use
       as an argument to call the jacobi function.
       Here are the matrix values stored in row-major
       order in a 1-D array. This array name will be used
       as an argument for the convert_matrix function.
     */
    float a_array[] = {
         5.0, -2.0, 
        -2.0,  2.0

    };

    float *d = vector(1, n);
    float *r = vector(1, n);
    float **v = matrix(1, n, 1, n);
    float **a = convert_matrix(&a_array[0], 1, n, 1, n);

    printf("Matrix:\n");
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            printf("  %+13.6e", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");


/*Jacobi*/

 

/*Fine jacobi*/

    //jacobi(a, n, d, v, &nrot);




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

   // printf("\n");
   // printf("tempo impiegato: %f ", MPI_Wtime() - w_time);

    return 0;
   
}


