#include <stdio.h>
#include <iostream>
#include <time.h>
#include "getmatrix.hpp"
#include "myalgo.hpp"
#include "printmatrix.hpp"
#include "mistakematrix.hpp"
#include "mpi.h"

void matrix_init(int n,  double * initmat); 
int get_matrix(int n, const char* filename, double* A);
void print_matrix(int n, const double* A); 
int my_algo(int n, double* A, double* ans);
double func(int k, int n, int i, int j); 
void calc_mistake(double* res, int size);
void compinterval(int start, int interval, double* A, double* B, double* AB, int N);



int main(int argc, char **argv)
{

    MPI_Status stat;
    int n = atoi(argv[1]);
    int res = 0;
    double *A = new double[n * n];
    double *temp = new double[n * n];
    double* ans = new double[n * n];
    int interval, remainder;
    int size, rank;
    double * temp2 = new double[n * n];
    double time1, time2; 
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {   
            temp2[i*n + j] = 0;
        }
    }
    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    interval = n/(size);
    remainder = n % (size);
    if (rank == 0)
    {
       // std::cerr << argc << "\n";
        if ((argc != 2) && (argc != 3)) {return -1;} 
        


        //std::cerr << n << "\n";

        if (argc == 2) {
            matrix_init(n, A);
            matrix_init(n, temp);
        } else {
            res = get_matrix(n, argv[2], A);
            res = get_matrix(n, argv[2], temp);
            if (res != 0)  {return -1;}
        }

       // print_matrix(n, A);

        clock_t tStart = clock();
        res = my_algo(n, A, ans);
        clock_t tEnd = clock();

        if (res != 0)  {return -1;}

        printf("Time taken: %.2fs\n", (double)(tEnd - tStart)/CLOCKS_PER_SEC);

        // print_matrix(n, ans);
        // printf("\n helooooo\n");
        // print_matrix(n, temp);

    
        time1 = MPI_Wtime();  
        MPI_Bcast(ans,n*n,MPI_DOUBLE,0,MPI_COMM_WORLD); // send broadcast
        printf("%d: Bcast complete\n",rank);

        // Send intervals of array A to worker processes
        for(int i=1;i<size;i++)
            MPI_Send(temp+(i*interval*n),interval*n,MPI_DOUBLE,i,i,MPI_COMM_WORLD);
      
        compinterval(0,interval, temp, ans, temp2, n);                // local work
        compinterval(size*interval, remainder, temp, ans, temp2, n);  // remainder
        printf("local work finished\n");
      //get results from workers
        for(int i=1;i<size;i++)
	        MPI_Recv(temp2+(i*interval*n),interval*n,MPI_DOUBLE,i,i, MPI_COMM_WORLD, &stat);
        // print_matrix(n, temp2);
        time2 = MPI_Wtime();  
        printf("approx %d-process time Tp: %f sec.\n",size,time2-time1);  
        delete[] A;
        delete[] ans;
        delete[] temp;

        calc_mistake(temp2, n);
    }else {
        MPI_Bcast(ans,n*n,MPI_DOUBLE,0,MPI_COMM_WORLD); // receive broadcast
      // synchronous receive

      MPI_Recv(temp+(rank*interval*n),interval*n,MPI_DOUBLE,0,rank,
	       MPI_COMM_WORLD,&stat);
      
      compinterval(rank*interval, interval, temp, ans, temp2, n);
      // send results back to root process, synchronous send
        printf("process %d finished, temp2 = %f\n", rank, temp2[rank*interval*n + rank*interval]);
      //print_matrix(interval, temp2+(rank*interval));

      MPI_Send(temp2+(rank*interval*n),interval*n,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
    }
    delete[] temp2;
 
    return 0;
}

void compinterval(int start, int interval, double* A, double* B, double* AB, int N)
{
  int i, j, k;
//   printf("start = %d, interval = %d\n", start, interval);
  for (i=start;i<start+interval;i++)
    for(j=0;j<N;j++)
    {
	    AB[i*N+j]=0;
        
	    for(k=0;k<N;k++)
        {
            // if(j == start)
            //     printf("%f*%f + ", A[i*N+k],B[k*N+j]);
            AB[i*N+j] += A[i*N +k] * B[k*N+j];
        } 
        // printf("\n");
       // printf("AB[%d][%d] = %f\n", i, j, AB[i*N+j]);
    }
    // printf("AB[%d][%d] = %f\n", start, start, AB[start*N + start]);
}