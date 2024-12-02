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

int main(int argc, char **argv)
{
    MPI_Status status;
    int n = atoi(argv[1]); 
    double* ans = new double[n * n];
    double* eigenval = new double[n];
    double precision = 0.000001;
    double * temp2 = new double[n * n];
    int rank, size;
    int res = 0;
    double *A = new double[n * n];
    double *temp = new double[n * n];
    
    int rowSize, lastRowSize;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
        lastRowSize = n%(n/(size-1));
        rowSize = (n - lastRowSize)/(size-1);
        printf("rowSize = %d, lastRowSize = %d, n, size = %d, %d\n", rowSize, lastRowSize, n, size);
        for(int to_thread = 1; to_thread < size; to_thread++)
        {
            MPI_Send(&A, n*n, MPI_DOUBLE, to_thread, 0, MPI_COMM_WORLD);
            MPI_Send(&ans, n*n, MPI_DOUBLE, to_thread, 0, MPI_COMM_WORLD);
            MPI_Send(&n, 1, MPI_INT, to_thread, 0, MPI_COMM_WORLD);
            MPI_Send(&rowSize, 1, MPI_INT, to_thread, 0, MPI_COMM_WORLD);
            MPI_Send(&lastRowSize, 1, MPI_INT, to_thread, 0, MPI_COMM_WORLD);
        }
    }else {
        MPI_Recv(&A, n*n, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    MPI_Recv(&ans, n*n, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	    MPI_Recv(&n, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&rowSize, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&lastRowSize, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }


    
    printf("rank = %d\n", rank);
    mpiMatrMult(A, ans, temp2, n, rank, rowSize);
    MPI_Finalize();
    calc_mistake(temp2, n);

    delete[] A;
    delete[] ans;
    delete[] temp;
    delete[] temp2;
 
    return 0;
}