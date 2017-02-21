
/*  ***********************************************************
 *  Serial program to solve a Laplaces equation 
 *  (d2/dx2)u + (d2/dy2)u  = 0  using Jacobi method
 * *********************************************************** */

#include <stdio.h>
#include <math.h>

#define N 200

void main()
{
  double u[N+2][N+2], unew[N+2][N+2];
  double error;
  int iter, i, j;

  for(i=1; i<=N ; i++)
    for(j=1; j<=N ; j++)
      {
	u[i][j] = 0.0;
	unew[i][j] = u[i][j];
      }


      
  /*
    Boundary conditions:
 
                   u=10
                 --------
                |        |
                |        |
            u=0 |        | u=0
                |        |
                |        |
                 --------
                   u=0
  */
  
  for(i=1; i<=N ; i++)
    {
      u[i][0] = 0.0;
      u[i][N+1] = 10.0;
    }
  for(j=1; j<=N ; j++)
    {
      u[0][j] = 0.0;
      u[N+1][j] = 0.0;
    }
  

  /* Main iterative loop */
  error = 1000.0;
  iter = 0;
  
  while ( error > 1.0e-6  && iter < 1000){
    
    for(i=1; i<=N ; i++)
      for(j=1; j<=N ; j++)
	unew[i][j] = 0.25*(  u[i+1][j] + u[i-1][j] 
			     + u[i][j+1] + u[i][j-1] );
    
    error = 0.0;
    for(i=1; i<=N ; i++)
      for(j=1; j<=N ; j++)
	{
	  error += fabs(unew[i][j]-u[i][j]);
	  u[i][j] = unew[i][j];
	}
         
   iter++;

    if( iter%100 == 0)  printf("Iteration %5d error = %e \n",iter,error);

  }
}



