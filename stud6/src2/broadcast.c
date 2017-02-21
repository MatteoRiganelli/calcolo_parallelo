#include <mpi.h>
#include <stdio.h>
int main(int argc, char *argv[]) {
           int id;
           int n;
           int var_to_share;
           MPI_Init(&argc, &argv);
           MPI_Comm_size(MPI_COMM_WORLD, &n);
           MPI_Comm_rank(MPI_COMM_WORLD, &id);
           if (id==0) {
                  var_to_share=5;
           }
           MPI_Bcast (&var_to_share,1,MPI_INT,0,MPI_COMM_WORLD);
           printf("Processo %d: Ho a disposizione il valore %d;\n",id, var_to_share);
           MPI_Finalize();
           return 0;
}