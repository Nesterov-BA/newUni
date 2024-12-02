#include <stdio.h>
#include <iostream>
#include <time.h>
#include "getmatrix.hpp"
#include "myalgo.hpp"
#include "printmatrix.hpp"
#include "mistakematrix.hpp"
#include "functions.hpp"
#include "matrspin.hpp"
#include "qrreflect.hpp"
#include "mpi.h"


void matrix_init(int n,  double * initmat); 
int get_matrix(int n, const char* filename, double* A);
void print_matrix(int n, const double* A); 
int my_algo(int n, double* A, double* ans);
double func(int k, int n, int i, int j); 
void calc_mistake(int n, double* A, double* ans, double* temp); 
void mpiMatrMult(double* A, double* B, double* res, int size, int rank, int rowSize);
void calc_mistake(double* res, int size);
void compinterval(int start, int interval, double* A, double* B, double* AB, int N);

int main(int argc, char **argv)
{
    MPI_Status stat;
    int n = atoi(argv[1]); 
    double* ans = new double[n * n];
    double* eigenval = new double[n];
    double precision = 0.000001;
    double * temp2 = new double[n * n];
    int rank, size;
    int res = 0;
    double *A = new double[n * n];
    double *temp = new double[n * n];
    int interval, remainder;
    
    int rowSize, lastRowSize;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    interval = n/(size);
    remainder = n % (size);
    if (rank == 0)
    {
        std::cerr << argc << "\n";
        if ((argc != 2) && (argc != 3)) {return -1;} 


        std::cerr << n << "\n";


        if (argc == 2) {
            matrix_init(n, A);
            matrix_init(n, temp);
        } else {
            res = get_matrix(n, argv[2], A);
            res = get_matrix(n, argv[2], temp);
            if (res != 0)  {return -1;}
        }

        //print_matrix(n, A);
        

        id(ans, n);

        clock_t tStart = clock();
        print_matrix(n, A);
        triang(A, n);
        printf("tridiag:\n");
        print_matrix(n, A);
        qrref(A, ans, eigenval, precision, n);
        clock_t tEnd = clock();


        printf("Time taken: %.2fs\n", (double)(tEnd - tStart)/CLOCKS_PER_SEC);

        print_matrix(n, A);
        printf("Тут наны\n");
        print_matrix(n, ans);


        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                temp2[i*n + j] = 0;
            }
        }
        if(n%size == 0)
        {

        }

        MPI_Bcast(ans,n*n,MPI_DOUBLE,0,MPI_COMM_WORLD); // send broadcast
        printf("%d: Bcast complete\n",rank);

        // Send intervals of array A to worker processes
        for(int i=1;i<size;i++)
            MPI_Send(A+(i*interval),interval*n,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
      
        compinterval(0,interval, A, ans, temp2, n);                // local work
        compinterval(size*interval, remainder, A, ans, temp2, n);  // remainder

      //get results from workers
        for(int i=1;i<size;i++)
	        MPI_Recv(temp2+(i*interval),interval*n,MPI_DOUBLE,i,i, MPI_COMM_WORLD, &stat);

        calc_mistake(temp2, n);
    }else {
        MPI_Bcast(ans,n*n,MPI_DOUBLE,0,MPI_COMM_WORLD); // receive broadcast
      // synchronous receive

        printf("process %d finished\n", rank);
      MPI_Recv(A+(rank*interval),interval*n,MPI_DOUBLE,0,rank,
	       MPI_COMM_WORLD,&stat);
      
      compinterval(rank*interval, interval, A, ans, temp2, n);
      // send results back to root process, synchronous send
      
      MPI_Send(temp2+(rank*interval),interval*n,MPI_DOUBLE,0,rank,
		MPI_COMM_WORLD);
    }


    
    // printf("rank = %d\n", rank);
    // mpiMatrMult(A, ans, temp2, n, rank, rowSize);
    // printf("hello, rank = %d\n", rank);
    MPI_Finalize();

    delete[] A;
    delete[] ans;
    delete[] temp;
    delete[] temp2;
 
    return 0;
}

void compinterval(int start, int interval, double* A, double* B, double* AB, int N)
{
  int i, j, k;
  for (i=start;i<start+interval;i++)
    for(j=0;j<N;j++)
      {
	AB[i*N+j]=0;
	for(k=0;k<N;k++) AB[i*N+j] += A[i*N +k] * B[k*N+j];
      }
}